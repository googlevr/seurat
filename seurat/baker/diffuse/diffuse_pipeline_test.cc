/*
Copyright 2017 Google Inc. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS-IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include "seurat/baker/diffuse/diffuse_pipeline.h"

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "absl/memory/memory.h"
#include "seurat/artifact/artifact.h"
#include "seurat/baker/diffuse/diffuse_baker.h"
#include "seurat/baker/framework/rasterizer.h"
#include "seurat/baker/framework/texture_sizer.h"
#include "seurat/base/color.h"
#include "seurat/geometry/binning_point_cloud_builder.h"
#include "seurat/geometry/point_cloud_builder.h"
#include "seurat/geometry/quad_mesh.h"
#include "seurat/image/color_processor.h"
#include "seurat/image/filter.h"
#include "seurat/image/image.h"
#include "seurat/image/image_util.h"
#include "seurat/image/nearly_square_atlaser.h"
#include "seurat/ingest/point_cloud_assembler.h"
#include "seurat/ingest/view_group_loader_test_utils.h"
#include "seurat/mesh/mesh_component.h"
#include "seurat/mesh/mesh_component_util.h"
#include "seurat/tiler/tiler.h"

namespace seurat {
namespace baker {
namespace diffuse {
namespace {

using artifact::Artifact;
using base::Color1f;
using base::Color3f;
using geometry::BinningPointCloudBuilder;
using geometry::PointCloudBuilder;
using geometry::QuadMesh;
using image::Atlaser;
using image::ColorProcessor;
using image::GammaToneMapper;
using image::Image4f;
using image::NearlySquareAtlaser;
using ingest::FakeViewGroupLoader;
using ingest::PointCloudAssembler;
using ingest::ViewGroupLoader;
using ion::math::Vector2i;
using ion::math::Vector3i;
using mesh::MeshComponent;
using mesh::MeshComponentUtil;
using tiler::Tiler;
using tiler::TilerFactory;

const int kThreadCount = 3;

std::unique_ptr<ViewGroupLoader> MakeTestingViewGroupLoader() {
  const int kNumViewGroups = 2;
  const Vector2i kImageSize(64, 64);
  std::array<Color3f, 6> face_colors = {{
      {2.0f, 0.0f, 0.0f},  //
      {0.0f, 1.0f, 0.0f},  //
      {0.0f, 0.0f, 1.0f},  //
      {1.0f, 1.0f, 0.0f},  //
      {0.0f, 1.0f, 1.0f},  //
      {0.0f, 0.0f, 0.0f}   //
  }};
  std::array<float, 6> face_depths = {{0.1f, 0.5f, 0.3f, 0.6f, 0.4f, 0.2f}};

  return std::unique_ptr<ViewGroupLoader>(new FakeViewGroupLoader(
      kNumViewGroups, kImageSize, face_colors, face_depths));
}

TEST(DiffusePipelineTest, Process) {
  const Vector3i kBinningResolution(2048, 2048, 1 << 16);
  const float kNearClip = 0.01f;
  std::unique_ptr<PointCloudBuilder> point_cloud_builder(
      new BinningPointCloudBuilder(kThreadCount, kBinningResolution,
                                   kNearClip));
  std::unique_ptr<ingest::PointCloudAssembler> merger(
      new PointCloudAssembler(kThreadCount, MakeTestingViewGroupLoader(),
                              std::move(point_cloud_builder)));

  const int kMaxTileCount = 64;
  TilerFactory::Parameters tiler_parameters;
  tiler_parameters.tile_count = kMaxTileCount;
  tiler_parameters.peak_overdraw_factor = 999.0f;
  tiler_parameters.peak_overdraw_samples = 0;
  tiler_parameters.overdraw_factor = 2.0f;
  tiler_parameters.thread_count = kThreadCount;
  tiler_parameters.min_subdivision_level = 1;
  tiler_parameters.max_subdivision_level = 2;
  tiler_parameters.fast = true;
  std::unique_ptr<tiler::Tiler> tiler =
      TilerFactory::CreateDefaultTiler(tiler_parameters);

  std::unique_ptr<FrameGenerator> frame_generator(new TilerFrameGenerator(
      std::move(merger), std::move(tiler),
      std::unique_ptr<baker::FrameSorter>(new baker::FrameSorter(kThreadCount)),
      TilerFrameGenerator::TextureParameterization::kObjectSpace));
  const float kSecondaryFrameThreshold = 0.015f;
  std::unique_ptr<baker::RayClassifier> ray_classifier(
      new baker::ProjectingRayClassifier(
          kThreadCount,
          baker::ProjectingRayClassifier::RenderingMode::kDrawOrder,
          kSecondaryFrameThreshold));
  std::unique_ptr<FrameRasterizer> rasterizer(new FrameRasterizer(
      kThreadCount, MakeTestingViewGroupLoader(), std::move(ray_classifier)));
  const float kPixelsPerDegree = 0.5f;
  std::unique_ptr<TextureSizer> texture_sizer(
      new ProjectedAreaTextureSizer(kPixelsPerDegree));
  const float kSpecularFilterSize = 0.10f;
  const float kSigmaPixel = 0.3f;
  const float kRadiusPixel = 1.5f;
  std::shared_ptr<image::GaussianFilter> pixel_filter =
      std::make_shared<image::GaussianFilter>(kSigmaPixel, kRadiusPixel);

  DiffusePipeline diffuse_pipeline(
      kThreadCount, std::move(frame_generator), std::move(texture_sizer),
      absl::make_unique<DiffuseBaker>(kThreadCount, std::move(rasterizer),
                                      kSpecularFilterSize, pixel_filter),
      absl::make_unique<GammaToneMapper>(1.0f));

  Artifact artifact;
  EXPECT_TRUE(diffuse_pipeline.Process(&artifact).ok());

  EXPECT_TRUE(artifact.quad_mesh);
  const QuadMesh& quad_mesh = *artifact.quad_mesh;

  // Simply assert that the final quad mesh has one quad for each partition.
  EXPECT_GE(kMaxTileCount, quad_mesh.quads.size());

  // Verifies the existence of non-black texels.
  bool non_black_texel_exists = false;
  for (const auto& texture : quad_mesh.textures) {
    if (non_black_texel_exists) break;
    for (const auto& pixel : texture) {
      if (pixel[0] != 0 || pixel[1] != 0 || pixel[2] != 0) {
        non_black_texel_exists = true;
        break;
      }
    }
  }
  EXPECT_TRUE(non_black_texel_exists);
}

}  // namespace
}  // namespace diffuse
}  // namespace baker
}  // namespace seurat
