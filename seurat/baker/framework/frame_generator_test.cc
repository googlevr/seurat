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

#include "seurat/baker/framework/frame_generator.h"

#include <array>
#include <memory>
#include <vector>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/color.h"
#include "seurat/base/file_system.h"
#include "seurat/geometry/binning_point_cloud_builder.h"
#include "seurat/geometry/point_cloud_builder.h"
#include "seurat/ingest/point_cloud_assembler.h"
#include "seurat/ingest/view_group_loader_test_utils.h"
#include "seurat/testing/test_flags.h"
#include "seurat/tiler/tiler.h"

namespace seurat {
namespace baker {
namespace {

using base::Color3f;
using base::FileSystem;
using geometry::BinningPointCloudBuilder;
using geometry::PointCloudBuilder;
using ingest::FakeViewGroupLoader;
using ingest::PointCloudAssembler;
using ingest::ViewGroupLoader;
using ion::math::Vector2i;
using ion::math::Vector3i;
using tiler::Tiler;
using tiler::TilerFactory;

TEST(FrameGeneratorTest, FrameFromPartitionAndDiskCache_Tiler) {
  const int kNumViewGroups = 2;
  const Vector2i kImageSize(64, 64);
  std::array<Color3f, 6> face_colors = {{
      //
      {1.0f, 0.0f, 0.0f},  //
      {0.0f, 1.0f, 0.0f},  //
      {0.0f, 0.0f, 1.0f},  //
      {1.0f, 1.0f, 0.0f},  //
      {0.0f, 1.0f, 1.0f},  //
      {0.0f, 0.0f, 0.0f}   //
  }};
  std::array<float, 6> face_depths = {{0.1f, 0.5f, 0.3f, 0.6f, 0.4f, 0.2f}};

  const int kThreadCount = 1;
  std::unique_ptr<ViewGroupLoader> fake_loader(new FakeViewGroupLoader(
      kNumViewGroups, kImageSize, face_colors, face_depths));
  const Vector3i kBinningResolution(2048, 2048, 1 << 16);
  const float kNearClip = 0.01f;
  std::unique_ptr<PointCloudBuilder> point_cloud_builder(
      new BinningPointCloudBuilder(kThreadCount, kBinningResolution,
                                   kNearClip));
  std::unique_ptr<ingest::PointCloudAssembler> merger(new PointCloudAssembler(
      kThreadCount, std::move(fake_loader), std::move(point_cloud_builder)));
  std::unique_ptr<FrameSorter> sorter(new FrameSorter(kThreadCount));

  TilerFactory::Parameters tiler_parameters;
  tiler_parameters.tile_count = 64;
  tiler_parameters.overdraw_factor = 2.0f;
  tiler_parameters.peak_overdraw_factor = 4.0f;
  tiler_parameters.thread_count = 2;
  tiler_parameters.min_subdivision_level = 2;
  tiler_parameters.max_subdivision_level = 3;
  tiler_parameters.fast = true;
  std::unique_ptr<Tiler> tiler =
      TilerFactory::CreateDefaultTiler(tiler_parameters);

  std::unique_ptr<FrameGenerator> tiler_frame_generator(new TilerFrameGenerator(
      std::move(merger), std::move(tiler), std::move(sorter),
      TilerFrameGenerator::TextureParameterization::kObjectSpace));

  auto files = std::make_shared<base::FileSystem>(testing::GetTestTmpdir());
  const char kCacheFilename[] =
      "FrameGeneratorTest_FrameFromPartitionAndDiskCache_Tiler";
  // Clear any existing file.
  ASSERT_TRUE(files->SetContents(kCacheFilename, "").ok());
  DiskCachingFrameGenerator frame_generator(kCacheFilename, files,
                                            std::move(tiler_frame_generator));

  std::vector<Frame> frames;
  EXPECT_TRUE(frame_generator.GenerateFrames(&frames).ok());
  EXPECT_GE(64, frames.size());

  std::vector<Frame> frames2;
  EXPECT_TRUE(frame_generator.GenerateFrames(&frames2).ok());
  EXPECT_EQ(frames, frames2);
}

}  // namespace
}  // namespace baker
}  // namespace seurat
