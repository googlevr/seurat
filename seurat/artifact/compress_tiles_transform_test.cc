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

#include "seurat/artifact/compress_tiles_transform.h"

#include <vector>

#include "ion/gfx/image.h"
#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "absl/types/span.h"
#include "seurat/artifact/artifact.h"
#include "seurat/artifact/artifact_test_util.h"
#include "seurat/base/array2d.h"
#include "seurat/base/array2d_util.h"
#include "seurat/base/color.h"
#include "seurat/compressor/rgba/rgba_compressor.h"
#include "seurat/image/image.h"
#include "seurat/image/image_test_utils.h"
#include "seurat/image/image_util.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace artifact {
namespace {

using base::Color4f;
using geometry::QuadMesh;
using image::FakeAtlaser;
using image::Image4f;
using ion::gfx::ImagePtr;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2i;

class FakeTileCompressor : public compressor::RgbaCompressor {
 public:
  FakeTileCompressor() {}

  // RgbaCompressor implementation
  void Compress(absl::Span<const image::Image4f> textures,
                absl::Span<image::Image4f> compressed_textures) const override {
    // Note this is a terrible downsampling algorithm!
    std::transform(textures.begin(), textures.end(),
                   compressed_textures.begin(), [](const Image4f& texture) {
                     Image4f resized(texture);
                     resized.Resize(resized.GetSize() / 2);
                     return resized;
                   });
  }
};

// Tests codec operates during the transform.
TEST(CompressTilesTransformTest, CompressedArtifactAltersTextureSize) {
  // Builds the transform function under test.
  CompressTilesTransform transform(
      std::unique_ptr<compressor::RgbaCompressor>(new FakeTileCompressor),
      std::make_shared<FakeAtlaser>());
  // Prepares input data.
  const Color4f kSlightRed{125.5f / 255.0f, 0.0f, 0.0f, 1.0f};

  Artifact quad_artifact;
  quad_artifact.quad_mesh =
      std::make_shared<QuadMesh>(MakeTwoQuadMesh(kSlightRed));
  const QuadMesh initial_quad_mesh = *quad_artifact.quad_mesh;

  // Exercises the transformation under test.
  ASSERT_TRUE(transform.Process(&quad_artifact).ok());
  const QuadMesh& final_quad_mesh = *quad_artifact.quad_mesh;

  // Validates final objects.
  for (auto const& quad : final_quad_mesh.quads) {
    EXPECT_EQ(final_quad_mesh.textures[quad.texture_index].GetSize(),
              initial_quad_mesh.textures[quad.texture_index].GetSize() / 2);
  }
}

// Tests bad input artifacts die.
TEST(CompressTilesTransformTest, BadInputArtifactChecks) {
  CompressTilesTransform transform(
      std::unique_ptr<compressor::RgbaCompressor>(new FakeTileCompressor),
      std::make_shared<FakeAtlaser>());
  // Prepares "bad" input data.
  Artifact bad_input_artifact;

  // Exercises the transformation under test.
  base::Status process_status;
  EXPECT_DEATH(
      process_status = transform.Process(&bad_input_artifact),
      "The input artifact does not support the required representation");
}

}  // namespace
}  // namespace artifact
}  // namespace seurat
