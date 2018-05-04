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

#include "seurat/base/array2d.h"

#include <numeric>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "absl/types/span.h"
#include "seurat/base/color.h"

namespace seurat {
namespace base {
namespace {

using base::Color4f;
using ion::math::Point2i;

// Basic construction tests for Array2D. Image types are aliases of Array2D.
// More unit tests are in image/image_test.cc.
TEST(Array2DTest, ConstructEmptyArray) {
  Array2D<Color4f> image;
  EXPECT_EQ(0, image.Width());
  EXPECT_EQ(0, image.Height());
}

// Copy test for Array2D.
TEST(Array2DTest, CopyTest) {
  const int kSize = 2;
  const Array2D<Color4f> copy_source(kSize, kSize * 2);

  // Check copy constructor
  Array2D<Color4f> copy_result(copy_source);
  EXPECT_EQ(copy_source, copy_result);

  // Check inequality.
  copy_result.Resize(0, 0);
  EXPECT_NE(copy_source, copy_result);

  // Check assignment operator.
  copy_result = copy_source;
  EXPECT_EQ(copy_source, copy_result);
}

// Move test for Array2D.
TEST(Array2DTest, MoveTest) {
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kGreen(0.0f, 1.0f, 0.0f, 1.0f);
  const int kSize = 2;
  Array2D<Color4f> move_source(kSize, kSize, kRed);
  EXPECT_EQ(move_source.Width(), kSize);
  EXPECT_EQ(move_source.Height(), kSize);
  EXPECT_NE(move_source.Data(), nullptr);

  // Check move construction.
  Array2D<Color4f> move_result(std::move(move_source));
  EXPECT_EQ(move_result.Width(), kSize);
  EXPECT_EQ(move_result.Height(), kSize);
  EXPECT_NE(move_result.Data(), nullptr);

  // Tests the object behaves safely after the move, including reinitialization
  // after the move.
  move_source.Resize(kSize, kSize, kGreen);  // NOLINT misc-use-after-move
  EXPECT_TRUE(move_source.IsInside(Point2i(0, 0)));
  EXPECT_EQ(move_source.Width(), kSize);
  EXPECT_EQ(move_source.Height(), kSize);

  // Check move assignment.
  move_result.Resize(0, 0);
  EXPECT_NE(move_source, move_result);
  move_result = std::move(move_source);
  EXPECT_EQ(move_result.Width(), kSize);
  EXPECT_EQ(move_result.Height(), kSize);

  Array2D<Color4f> copy_result(move_result);
  EXPECT_EQ(copy_result, move_result);
}

TEST(Array2DTest, ConstructArrayWithDefaultValues) {
  const Color3f kZero(0.0f, 0.0f, 0.0f);
  Array2D<Color3f> image(2, 3);
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(3, image.Height());
  for (int y = 0; y < image.Height(); y++) {
    for (int x = 0; x < image.Width(); x++) {
      EXPECT_EQ(kZero, image.At(x, y));
    }
  }
}

TEST(Array2DTest, ConstructArrayWithGivenValues) {
  const Color3f kTestValue(0.0f, 0.2f, 0.4f);
  Array2D<Color3f> image(2, 3, kTestValue);
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(3, image.Height());
  for (int y = 0; y < image.Height(); y++) {
    for (int x = 0; x < image.Width(); x++) {
      EXPECT_EQ(kTestValue, image.At(x, y));
    }
  }
}

TEST(Array2DTest, ConstructBoolArrayWithGivenValues) {
  const bool kTestValue(true);
  Array2D<bool> bitmap(2, 3, kTestValue);
  EXPECT_EQ(2, bitmap.Width());
  EXPECT_EQ(3, bitmap.Height());
  for (int y = 0; y < bitmap.Height(); y++) {
    for (int x = 0; x < bitmap.Width(); x++) {
      EXPECT_EQ(kTestValue, bitmap.At(x, y));
    }
  }
}

TEST(Array2DTest, BeginEnd) {
  Array2D<int> image(3, 2);
  std::iota(image.begin(), image.end(), 0);
  EXPECT_EQ(0, image.At(0, 0));
  EXPECT_EQ(5, image.At(2, 1));
}

TEST(Array2DTest, Span) {
  Array2D<int> image(3, 2);
  absl::Span<int> image_span = absl::MakeSpan(image);
  std::iota(image_span.begin(), image_span.end(), 0);
  EXPECT_EQ(0, image.At(0, 0));
  EXPECT_EQ(5, image.At(2, 1));
}

}  // namespace
}  // namespace base
}  // namespace seurat
