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

#include "seurat/image/image.h"

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/color.h"

using ion::math::Point2i;
using seurat::base::Color4f;

namespace seurat {
namespace image {
namespace {

TEST(ImageTest, AtTest) {
  Image4f image(2, 2);

  const Color4f kZero(0.0f, 0.0f, 0.0f, 0.0f);
  const Color4f kTest(1.0f, 4.0f, 9.0f, 25.0f);
  image.Fill(kZero);
  EXPECT_EQ(kZero, image.At({0, 1}));
  image.At(0, 1) = kTest;
  EXPECT_EQ(kTest, image.At({0, 1}));

  // Test different coordinates are accessible but not modified by prior
  // assignment.
  EXPECT_EQ(kZero, image.At({0, 0}));
  EXPECT_EQ(kZero, image.At({1, 0}));
  EXPECT_EQ(kZero, image.At({1, 1}));

#if defined(ION_BLAZE) && defined(ION_PLATFORM_LINUX)
  // EXPECT_DEBUG_DEATH does not work properly on Windows.
  EXPECT_DEBUG_DEATH(image.At({image.Width(), 0}), "IsInside");
#endif
}

TEST(ImageTest, EqualityTest) {
  Image4f image1(16, 8);
  Image4f image2(8, 16);
  Image4f image3(16, 8);
  Image4f image3_copy(16, 8);

  image1.Fill({1.0f, 2.0f, 3.0f, 4.0f});
  image2.Fill({1.0f, 2.0f, 3.0f, 4.0f});
  image3.Fill({4.0f, 3.0f, 2.0f, 1.0f});
  image3_copy.Fill({4.0f, 3.0f, 2.0f, 1.0f});

  EXPECT_FALSE(image1 == image2);
  EXPECT_FALSE(image3 == image1);
  EXPECT_TRUE(image1 == image1);
  EXPECT_TRUE(image2 == image2);
  EXPECT_TRUE(image3 == image3);
  EXPECT_TRUE(image3_copy == image3);
}

TEST(ImageTest, BoundsTest) {
  Image4f image(16, 8);

  const Point2i kInsideCorners[] = {{0, 0}, {15, 0}, {15, 7}, {0, 7}};
  for (auto p : kInsideCorners) {
    EXPECT_TRUE(image.IsInside(p)) << p;
  }

  const Point2i kOutsideCorners[] = {{-1, 0}, {0, -1},  {-1, -1},  // L, B, LL
                                     {16, 0}, {16, -1}, {15, -1},  // R, B, LR
                                     {16, 7}, {15, 8},  {16, 8},   // R, T, UR
                                     {-1, 7}, {0, 8},   {-1, 8}};  // L, T, UL
  for (auto p : kOutsideCorners) {
    EXPECT_FALSE(image.IsInside(p)) << p;
  }
}

TEST(ImageTest, ImageBool) {
  Image<bool> image(2, 2);

  image.Fill(false);
  EXPECT_EQ(false, image.At(0, 1));
  EXPECT_EQ(false, image.At(1, 1));
  EXPECT_EQ(false, image.At(1, 0));
  EXPECT_EQ(false, image.At(1, 1));

  image.At(0, 1) = true;
  EXPECT_EQ(true, image.At(0, 1));

  image.At(0, 1) = false;
  EXPECT_EQ(false, image.At(0, 1));
}

TEST(ImageTest, DataTest) {
  // Data should be in row-major order without padding.
  Image4f image_4f(2, 2);
  image_4f.Fill(Color4f::Zero());

  image_4f.At(0, 1) = Color4f(1.0f, 2.0f, 3.0f, 4.0f);
  EXPECT_EQ(Color4f(1.0f, 2.0f, 3.0f, 4.0f), image_4f.Data()[2]);

  image_4f.Data()[2] = Color4f(42.0f, 42.0f, 42.0f, 42.0f);
  EXPECT_EQ(Color4f(42.0f, 42.0f, 42.0f, 42.0f), image_4f.At(0, 1));
}

}  // namespace
}  // namespace image
}  // namespace seurat
