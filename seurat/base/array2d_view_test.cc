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

#include "seurat/base/array2d_view.h"

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/array2d.h"
#include "seurat/base/array2d_util.h"
#include "seurat/base/color.h"

using ion::math::Point2i;
using seurat::base::Color4f;

namespace seurat {
namespace base {
namespace {

// Basic construction tests for Array2DView.
TEST(Array2DViewTest, ConstructEmptyArrayView) {
  Array2DView<Color4f> view;
  EXPECT_EQ(0, view.Width());
  EXPECT_EQ(0, view.Height());
}

TEST(Array2DViewTest, ConstructFromArray2D) {
  const Color3f kRed(1.0f, 0.0f, 0.0f);
  Array2D<Color3f> arr(2, 3, kRed);
  Array2DView<Color3f> view(&arr);
  EXPECT_EQ(2, view.Width());
  EXPECT_EQ(3, view.Height());
  for (int y = 0; y < view.Height(); ++y) {
    for (int x = 0; x < view.Width(); ++x) {
      EXPECT_EQ(kRed, view.At(x, y));
    }
  }
}

Array2D<int> MakeTestArray() {
  Array2D<int> arr(5, 4);
  for (int y = 0; y < arr.Height(); ++y) {
    for (int x = 0; x < arr.Width(); ++x) {
      arr.At(x, y) = x + y * 10;
    }
  }

  return arr;
}

TEST(Array2DViewTest, ConstructFromArray2D_IdentityCrop) {
  Array2D<int> arr = MakeTestArray();
  int crop_width = arr.Width();
  int crop_height = arr.Height();
  Array2DView<int> view(&arr, 0, 0, crop_width, crop_height);
  EXPECT_EQ(crop_width, view.Width());
  EXPECT_EQ(crop_height, view.Height());
  for (int y = 0; y < view.Height(); ++y) {
    for (int x = 0; x < view.Width(); ++x) {
      EXPECT_EQ(x + y * 10, view.At(x, y));
    }
  }
}

TEST(Array2DViewTest, ConstructFromArray2D_InvalidCrop) {
// Expect DCHECK for dbg and fastbuild, on platforms where GTest supports it.
#if GTEST_HAS_DEATH_TEST
  Array2D<int> arr = MakeTestArray();
  EXPECT_DEATH(Array2DView<int>(&arr, arr.Width(), arr.Height(), 1, 1),
               "IsInside");

  // x & y are out of bounds
  EXPECT_DEATH(
      Array2DView<int>(&arr, arr.Width() + 10, arr.Height() + 10, 1, 1),
      "IsInside");

  EXPECT_DEATH(Array2DView<int>(&arr, 0, 0, -1, 3), "Invalid width");
  EXPECT_DEATH(Array2DView<int>(&arr, 0, 0, 100, 3), "Invalid crop");
  EXPECT_DEATH(Array2DView<int>(&arr, 0, 0, 1, -2), "Invalid height");
  EXPECT_DEATH(Array2DView<int>(&arr, 0, 0, 1, 100), "Invalid crop");
#endif
}

TEST(Array2DViewTest, ConstructFromArray2D_CropBottomRight) {
  Array2D<int> arr = MakeTestArray();
  Array2DView<int> view(&arr, 0, 0, 2, 3);
  EXPECT_EQ(2, view.Width());
  EXPECT_EQ(3, view.Height());
  for (int y = 0; y < view.Height(); ++y) {
    for (int x = 0; x < view.Width(); ++x) {
      EXPECT_EQ(x + y * 10, view.At(x, y));
    }
  }
}

TEST(Array2DViewTest, ConstructFromArray2D_CropTopLeft) {
  Array2D<int> arr = MakeTestArray();
  Array2DView<int> view(&arr, 2, 3, arr.Width() - 2, arr.Height() - 3);
  for (int y = 0; y < view.Height(); ++y) {
    for (int x = 0; x < view.Width(); ++x) {
      EXPECT_EQ((x + 2) + (y + 3) * 10, view.At(x, y));
    }
  }
}

TEST(Array2DViewTest, ConstructFromArray2D_DoubleCrop) {
  Array2D<int> arr = MakeTestArray();
  Array2DView<int> view1(&arr, 3, 2, 2, 2);
  EXPECT_EQ(2, view1.Width());
  EXPECT_EQ(2, view1.Height());
  Array2DView<int> view2(view1, 1, 1, 1, 1);
  EXPECT_EQ(1, view2.Width());
  EXPECT_EQ(1, view2.Height());
  EXPECT_EQ(&arr.At(3 + 1, 2 + 1), &view2.At(0, 0));
}

TEST(Array2DViewTest, AtTest) {
  Array2D<Color4f> arr(2, 2);
  Array2DView<Color4f> view(&arr);

  const Color4f kZero(0.0f, 0.0f, 0.0f, 0.0f);
  const Color4f kTest(1.0f, 4.0f, 9.0f, 25.0f);
  arr.Fill(kZero);
  EXPECT_EQ(kZero, view.At({0, 1}));
  arr.At(0, 1) = kTest;
  EXPECT_EQ(kTest, view.At({0, 1}));

  // Test different coordinates are accessible but not modified by prior
  // assignment.
  EXPECT_EQ(kZero, view.At({0, 0}));
  EXPECT_EQ(kZero, view.At({1, 0}));
  EXPECT_EQ(kZero, view.At({1, 1}));

// Expect DCHECK for dbg and fastbuild, on platforms where GTest supports it.
#if GTEST_HAS_DEATH_TEST
  EXPECT_DEBUG_DEATH(view.At({view.Width(), 0}), "IsInside");
#endif
}

TEST(Array2DViewTest, EqualityTest) {
  Array2D<Color4f> arr1(16, 8);
  Array2D<Color4f> arr2(8, 16);
  Array2D<Color4f> arr3(16, 8);

  Array2DView<Color4f> view1(&arr1);
  Array2DView<Color4f> view2(&arr2);
  Array2DView<Color4f> view3(&arr3);

  arr1.Fill({1.0f, 2.0f, 3.0f, 4.0f});
  arr2.Fill({1.0f, 2.0f, 3.0f, 4.0f});
  arr3.Fill({4.0f, 3.0f, 2.0f, 1.0f});

  EXPECT_NE(arr1, arr2);
  EXPECT_FALSE(arr1 == arr2);
  EXPECT_FALSE(arr3 == arr1);
  EXPECT_EQ(arr1, arr1);
  EXPECT_TRUE(arr1 == arr1);
  EXPECT_TRUE(arr2 == arr2);
  EXPECT_TRUE(arr3 == arr3);
  EXPECT_NE(view1, view2);
  EXPECT_FALSE(view1 == view2);
  EXPECT_FALSE(view3 == view1);
  EXPECT_EQ(view1, view1);
  EXPECT_TRUE(view1 == view1);
  EXPECT_TRUE(view2 == view2);
  EXPECT_TRUE(view3 == view3);
}

TEST(Array2DViewTest, EqualityTest_Cropped) {
  // Test equality of cropped subregions of two 2D arrays which are different
  // outside the tested subregions.
  Array2D<int> arr1(4, 5);
  arr1.Fill(1);
  for (int y = 1; y < 3; ++y) {
    for (int x = 1; x < 3; ++x) {
      arr1.At(x, y) = 2;
    }
  }

  Array2D<int> arr2(5, 7);
  arr2.Fill(3);
  for (int y = 2; y < 4; ++y) {
    for (int x = 2; x < 4; ++x) {
      arr2.At(x, y) = 2;
    }
  }

  Array2DView<int> view1(&arr1, 1, 1, 2, 2);
  Array2DView<int> view2(&arr2, 2, 2, 2, 2);

  EXPECT_TRUE(view1 == view2);
}

TEST(Array2DViewTest, FillTest_Cropped) {
  Array2D<int> arr(5, 4);
  arr.Fill(1);

  MutableArray2DView<int> view(&arr, 2, 1, 2, 2);
  view.Fill(2);

  // Expected pattern:
  //
  // 11111
  // 11221
  // 11221
  // 11111
  for (int y = 0; y < 4; ++y) {
    for (int x = 0; x < 5; ++x) {
      if (x >= 2 && x < 4 && y >= 1 && y < 3) {
        EXPECT_EQ(2, arr.At(x, y)) << "x = " << x << " y = " << y;
      } else {
        EXPECT_EQ(1, arr.At(x, y)) << "x = " << x << " y = " << y;
      }
    }
  }
}

TEST(Array2DViewTest, BoundsTest) {
  Array2D<float> arr(16, 8);
  Array2DView<float> view(&arr);

  const Point2i kInsideCorners[] = {{0, 0}, {15, 0}, {15, 7}, {0, 7}};
  for (auto p : kInsideCorners) {
    EXPECT_TRUE(view.IsInside(p)) << p;
  }

  const Point2i kOutsideCorners[] = {{-1, 0}, {0, -1},  {-1, -1},  // L, B, LL
                                     {16, 0}, {16, -1}, {15, -1},  // R, B, LR
                                     {16, 7}, {15, 8},  {16, 8},   // R, T, UR
                                     {-1, 7}, {0, 8},   {-1, 8}};  // L, T, UL
  for (auto p : kOutsideCorners) {
    EXPECT_FALSE(view.IsInside(p)) << p;
  }
}

TEST(Array2DViewTest, MutableToImmutable) {
  Array2D<float> arr(16, 8);
  MutableArray2DView<float> mutable_view(&arr);

  base::SpatialForEachArrayEntry(
      mutable_view, [&mutable_view](const Point2i& p, const float& e) {
        mutable_view.At(p) = p[0] + mutable_view.Width() * p[1];
      });

  Array2DView<float> const_view(mutable_view);

  EXPECT_EQ(mutable_view.Width(), const_view.Width());
  EXPECT_EQ(mutable_view.Height(), const_view.Height());
  EXPECT_EQ(mutable_view.Stride(), const_view.Stride());
  base::SpatialForEachArrayEntry(
      mutable_view, [&const_view](const Point2i& p, const float& e) {
        EXPECT_EQ(const_view.At(p), e) << p;
      });
}

}  // namespace
}  // namespace base
}  // namespace seurat
