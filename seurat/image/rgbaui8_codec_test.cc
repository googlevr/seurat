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

#include "seurat/image/rgbaui8_codec.h"

#include "gtest/gtest.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace image {
namespace {

using base::Color4f;

// Tests RgbaUi8Codec is stable on empty images.
TEST(RgbaUi8CodecTest, EmptyImages) {
  RgbaUi8Codec codec;
  Image4f initial_image;
  ion::gfx::ImagePtr ion_image = codec.Compress(initial_image);
  Image4f final_image = codec.Decompress(ion_image);

  EXPECT_EQ(ion::gfx::Image::kRgba8888, codec.GetFormat());
  EXPECT_EQ(initial_image, final_image);
}

// Verifies RgbaUi8Codec is approximately an identity codec; quantization is
// expected.
TEST(RgbaUi8CodecTest, QuantizedIdentity) {
  const float kUi8Quantization = 0.5f / 255.0f;
  RgbaUi8Codec codec;
  Image4f initial_image(3, 2);
  initial_image.At(0, 0) = Color4f::Zero();
  initial_image.At(1, 0).Set(1.0f, 0.0f, 0.0f, 1.0f);
  initial_image.At(2, 0).Set(0.0f, 1.0f, 0.0f, 1.0f);
  initial_image.At(0, 1).Set(0.0f, 0.0f, 1.0f, 1.0f);
  initial_image.At(1, 1).Set(0.5f, 0.5f, 0.5f, 1.0f);
  initial_image.At(2, 1).Set(0.75f, 0.25f, 0.0f, 1.0f);

  ion::gfx::ImagePtr ion_image = codec.Compress(initial_image);

  Image4f final_image = codec.Decompress(ion_image);

  EXPECT_EQ(ion::gfx::Image::kRgba8888, codec.GetFormat());

  EXPECT_EQ(initial_image.At(0, 0), final_image.At(0, 0));
  EXPECT_EQ(initial_image.At(1, 0), final_image.At(1, 0));
  EXPECT_EQ(initial_image.At(2, 0), final_image.At(2, 0));

  EXPECT_EQ(initial_image.At(0, 1), final_image.At(0, 1));
  EXPECT_VECTOR_NEAR(initial_image.At(1, 0), final_image.At(1, 0),
                     kUi8Quantization);
  EXPECT_VECTOR_NEAR(initial_image.At(2, 0), final_image.At(2, 0),
                     kUi8Quantization);
}

}  // namespace
}  // namespace image
}  // namespace seurat
