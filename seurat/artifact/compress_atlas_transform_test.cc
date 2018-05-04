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

#include "seurat/artifact/compress_atlas_transform.h"

#include <vector>

#include "ion/gfx/image.h"
#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/artifact/artifact.h"
#include "seurat/artifact/artifact_test_util.h"
#include "seurat/artifact/atlas_mesh_transform.h"
#include "seurat/base/array2d.h"
#include "seurat/base/array2d_util.h"
#include "seurat/base/color.h"
#include "seurat/image/codec.h"
#include "seurat/image/image.h"
#include "seurat/image/image_test_utils.h"
#include "seurat/image/image_util.h"
#include "seurat/image/nearly_square_atlaser.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace artifact {
namespace {

using base::Color4f;
using image::Image4f;
using ion::gfx::ImagePtr;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2i;

class FakeCodec : public image::Codec {
 public:
  FakeCodec() {}

  // Codec implementation
  ion::gfx::Image::Format GetFormat() const override {
    return ion::gfx::Image::kRgba8;
  }
  Vector2i GetBlockSize() const override;
  ImagePtr Compress(const image::Image4f& image) const override;
  Image4f Decompress(const ion::gfx::ImagePtr& image) const override;
};

Vector2i FakeCodec::GetBlockSize() const { return Vector2i(1, 1); }

ImagePtr FakeCodec::Compress(const image::Image4f& image) const {
  return image::ConvertSeuratImageToIonImage(image);
}

Image4f FakeCodec::Decompress(const ion::gfx::ImagePtr& image) const {
  image::Image4ui8 rgbaui8_image =
      image::ConvertIonImageToSeuratImage<image::Image4ui8>(image);
  Image4f decompressed_image(rgbaui8_image.GetSize());
  base::TransformArray(rgbaui8_image, &decompressed_image,
                       [](const base::Color4ui8& c) { return c.AsColorF(); });
  return decompressed_image;
}

// Tests codec operates during the transform.
TEST(CompressAtlasTransformTest, BadArtifactTransformInput) {
  // Builds the transform function under test.
  CompressAtlasTransform transform(
      std::unique_ptr<image::Codec>(new FakeCodec));

  // Prepares input data.
  const Color4f kSlightRed{125.5f / 255.0f, 0.0f, 0.0f, 1.0f};
  Artifact artifact;
  artifact.quad_mesh =
      std::make_shared<geometry::QuadMesh>(MakeSingleQuadMesh(kSlightRed));

  AtlasMeshTransform atlas_mesh_transform(
      std::make_shared<image::NearlySquareAtlaser>());
  ASSERT_TRUE(atlas_mesh_transform.Process(&artifact).ok());
  ASSERT_TRUE(artifact.texture);
  const Image4f initial_texture = *artifact.texture;

  // Exercises the transformation under test.
  Artifact final_artifact = artifact;
  ASSERT_TRUE(transform.Process(&final_artifact).ok());

  ASSERT_TRUE(final_artifact.texture);
  const Image4f& final_texture = *final_artifact.texture;

  // Validates resulting color atlas is the compressed version of the input
  // atlas.
  base::SpatialForEachArrayEntry(
      final_texture,
      [&initial_texture](const Point2i& p, const Color4f& final_color) {
        const float kColorIntensityEpsilon = 0.5f / 255.0f;
        EXPECT_VECTOR_NEAR(final_color, initial_texture.At(p),
                           kColorIntensityEpsilon)
            << p;
      });
}

// Tests failing transforms die.
TEST(CompressAtlasTransformTest, IncompatibleInputArtifactDies) {
  // Builds the transform function under test.
  CompressAtlasTransform transform(
      std::unique_ptr<image::Codec>(new FakeCodec));
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
