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

#include "seurat/image/inpainting.h"

#include <array>
#include <random>

#include "gtest/gtest.h"
#include "seurat/base/array2d.h"
#include "seurat/base/color.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace image {
namespace {

using base::Array2D;
using base::Color4f;
using ion::math::Point2i;

constexpr float kEpsilon = 1.0e-6f;

TEST(InpaintingTest, TinyImages) {
  // Test with 0x0 and 1x1 images.
  Image4f image;
  Array2D<bool> mask;

  image.Resize(0, 0);
  mask.Resize(0, 0);
  // This should not crash.
  InpaintSmooth(mask, &image);

  image.Resize(1, 1);
  mask.Resize(1, 1);
  mask.At(0, 0) = false;
  image.At(0, 0) = {1.0f, 0.75f, 0.0f, 0.5f};
  InpaintSmooth(mask, &image);
  EXPECT_EQ(Color4f(1.0f, 0.75f, 0.0f, 0.5f), image.At(0, 0));

  image.Resize(1, 1);
  mask.Resize(1, 1);
  mask.At(0, 0) = true;
  image.At(0, 0) = {1.0f, 0.75f, 0.0f, 0.5f};
  InpaintSmooth(mask, &image);
  EXPECT_EQ(Color4f(1.0f, 0.75f, 0.0f, 0.5f), image.At(0, 0));
}

Image4f MakeTesting2x2Image() {
  Image4f image(2, 2);
  image.At(0, 0) = {1.0f, 0.0f, 0.0f, 0.0f};
  image.At(1, 0) = {0.0f, 1.0f, 0.0f, 0.0f};
  image.At(0, 1) = {0.0f, 0.0f, 1.0f, 0.0f};
  image.At(1, 1) = {0.0f, 0.0f, 0.0f, 1.0f};
  return image;
}

TEST(InpaintingTest, Inpaint2x2) {
  const Image4f kOriginalImage = MakeTesting2x2Image();
  {
    Image4f image = kOriginalImage;
    Array2D<bool> mask(2, 2);

    mask.Fill(false);
    mask.At(0, 0) = true;

    InpaintSmooth(mask, &image);

    const Color4f expected_color =
        (kOriginalImage.At(1, 0) + kOriginalImage.At(0, 1)) * 0.5f;

    EXPECT_VECTOR_NEAR(expected_color, image.At(0, 0), kEpsilon);
    EXPECT_VECTOR_NEAR(kOriginalImage.At(1, 0), image.At(1, 0), kEpsilon);
    EXPECT_VECTOR_NEAR(kOriginalImage.At(0, 1), image.At(0, 1), kEpsilon);
    EXPECT_VECTOR_NEAR(kOriginalImage.At(1, 1), image.At(1, 1), kEpsilon);
  }

  {
    Image4f image = kOriginalImage;
    Array2D<bool> mask(2, 2);

    mask.Fill(false);
    mask.At(0, 1) = true;

    InpaintSmooth(mask, &image);

    const Color4f expected_color =
        (kOriginalImage.At(0, 0) + kOriginalImage.At(1, 1)) * 0.5f;

    EXPECT_VECTOR_NEAR(kOriginalImage.At(0, 0), image.At(0, 0), kEpsilon);
    EXPECT_VECTOR_NEAR(kOriginalImage.At(1, 0), image.At(1, 0), kEpsilon);
    EXPECT_VECTOR_NEAR(expected_color, image.At(0, 1), kEpsilon);
    EXPECT_VECTOR_NEAR(kOriginalImage.At(1, 1), image.At(1, 1), kEpsilon);
  }
}

TEST(InpaintingTest, InpaintLargeImageWithIsolatedPixels) {
  // Large enough to perform several recursive iterations.
  Image4f large_image(13, 16);

  std::mt19937 random;
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  for (auto& color : large_image) {
    color = {dist(random), dist(random), dist(random), dist(random)};
  }

  Image4f inpainted_image = large_image;
  Array2D<bool> mask(large_image.GetSize());
  std::array<Point2i, 3> masked_pixels = {{{4, 7}, {8, 8}, {12, 0}}};
  for (const auto& coord : masked_pixels) {
    mask.At(coord) = true;
  }

  InpaintSmooth(mask, &inpainted_image);

  // Test that the mask was respected and all unmasked pixels were left
  // unmodified.
  for (int y = 0; y < large_image.Height(); ++y) {
    for (int x = 0; x < large_image.Width(); ++x) {
      if (!mask.At(x, y)) {
        EXPECT_EQ(large_image.At(x, y), inpainted_image.At(x, y));
      }
    }
  }
}

TEST(InpaintingTest, InpaintLargeImageWithLinearGradient) {
  // Large enough to perform several recursive iterations.
  Image4f large_image(13, 16);

  std::mt19937 random;
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  for (auto& color : large_image) {
    color = {dist(random), dist(random), dist(random), dist(random)};
  }

  const Color4f kTopColor(1.0f, 0.0f, 1.0f, 0.0f);
  const Color4f kBottomColor(0.0f, 1.0f, 1.0f, 1.0f);

  Array2D<bool> mask(large_image.GetSize());
  mask.Fill(true);
  // Set the top & bottom rows to a fixed color & inpaint everything between.
  for (int x = 0; x < large_image.Width(); ++x) {
    large_image.At(x, 0) = kBottomColor;
    large_image.At(x, large_image.Height() - 1) = kTopColor;
    mask.At(x, 0) = false;
    mask.At(x, large_image.Height() - 1) = false;
  }

  Image4f inpainted_image = large_image;
  InpaintSmooth(mask, &inpainted_image);

  // All rows should be the same.
  for (int y = 0; y < large_image.Height(); ++y) {
    for (int x = 0; x < large_image.Width(); ++x) {
      EXPECT_VECTOR_NEAR(inpainted_image.At(0, y), inpainted_image.At(x, y),
                         1e-3f);
    }
  }

  // Verify that a vertical gradient is generated.
  for (int y = 1; y < large_image.Height(); ++y) {
    for (int x = 0; x < large_image.Width(); ++x) {
      Color4f above_color = inpainted_image.At(x, y);
      Color4f below_color = inpainted_image.At(x, y - 1);
      EXPECT_GT(above_color[0], below_color[0]);
      EXPECT_LT(above_color[1], below_color[1]);
      EXPECT_EQ(above_color[2], below_color[2]);
      EXPECT_LT(above_color[3], below_color[3]);
    }
  }
}

TEST(InpaintingTest, Inpaint_NaN_LargeImageWithLinearGradient) {
  // Large enough to perform several recursive iterations.
  Image4f large_image(13, 16);

  // Start with a NaN image.
  Color4f nan_color;
  nan_color.Fill(std::numeric_limits<float>::quiet_NaN());
  large_image.Fill(nan_color);

  const Color4f kTopColor(1.0f, 0.0f, 1.0f, 0.0f);
  const Color4f kBottomColor(0.0f, 1.0f, 1.0f, 1.0f);

  Array2D<bool> mask(large_image.GetSize());
  mask.Fill(true);
  // Set the top & bottom rows to a fixed color & inpaint everything between.
  for (int x = 0; x < large_image.Width(); ++x) {
    large_image.At(x, 0) = kBottomColor;
    large_image.At(x, large_image.Height() - 1) = kTopColor;
    mask.At(x, 0) = false;
    mask.At(x, large_image.Height() - 1) = false;
  }

  Image4f inpainted_image = large_image;
  InpaintSmooth(mask, &inpainted_image);

  // All rows should be the same.
  for (int y = 0; y < large_image.Height(); ++y) {
    for (int x = 0; x < large_image.Width(); ++x) {
      EXPECT_VECTOR_NEAR(inpainted_image.At(0, y), inpainted_image.At(x, y),
                         1e-3f);
    }
  }

  // Verify that a vertical gradient is generated.
  for (int y = 1; y < large_image.Height(); ++y) {
    for (int x = 0; x < large_image.Width(); ++x) {
      Color4f above_color = inpainted_image.At(x, y);
      Color4f below_color = inpainted_image.At(x, y - 1);
      EXPECT_GT(above_color[0], below_color[0]);
      EXPECT_LT(above_color[1], below_color[1]);
      EXPECT_EQ(above_color[2], below_color[2]);
      EXPECT_LT(above_color[3], below_color[3]);
    }
  }
}

}  // namespace
}  // namespace image
}  // namespace seurat
