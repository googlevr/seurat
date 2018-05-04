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

#include "seurat/base/array2d_util.h"

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/array2d.h"
#include "seurat/base/array2d_view.h"
#include "seurat/base/color.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace base {
namespace {

using ion::math::Point2i;
using ion::math::Range1f;
using ion::math::Range2i;
using ion::math::Vector2i;

// Note these names emulate expected use in outer layers of libraries, e.g.
// arrays of color vectors. This is primarily to ease porting the test to base.

// Verify operations from source into destination produce size mismatch DCHECKs.
TEST(Array2dUtilTest, ArrayBoundsChecks) {
  Array2D<Color4f> image1(4, 3);
  Array2D<Color4f> dest;

  // Two passes: one with an empty destination, and a second with an non-empty
  // but still mismatched destination.
  for (int run_twice = 2; run_twice; --run_twice) {
    EXPECT_DEATH(CopyArray(image1, &dest, Vector2i(0, 0)), "Width\\(\\)");

    dest.Resize(image1.Width() - 1, image1.Height());
  }
}

TEST(Array2dUtilTest, CopyArray) {
  Array2D<Color4f> image1(4, 3);
  for (int y = 0; y < 3; ++y) {
    for (int x = 0; x < 4; ++x) {
      image1.At(x, y) = {0.0f, 0.5f, 1.0f, 0.5f};
    }
  }

  Array2D<Color4f> image2(7, 8);
  CopyArray(image1, &image2, {1, 1});

  for (int y = 0; y < image2.Height(); ++y) {
    for (int x = 0; x < image2.Width(); ++x) {
      if (y >= 1 && y < 4 && x >= 1 && x < 5) {
        EXPECT_EQ(0.0f, image2.At(x, y)[0]);
        EXPECT_EQ(0.5f, image2.At(x, y)[1]);
        EXPECT_EQ(1.0f, image2.At(x, y)[2]);
        EXPECT_EQ(0.5f, image2.At(x, y)[3]);
      } else {
        for (int c = 0; c < 4; ++c) {
          EXPECT_EQ(0.0f, image2.At(x, y)[c]);
        }
      }
    }
  }
}

// Clips an array by a range that has an empty intersection with the extent of
// the array.
TEST(Array2dUtilTest, ClipByDisjointRange) {
  const Vector2i kSize(16, 8);
  const Point2i kClipRangeMinPoint = Point2i(1, 2) + kSize;

  Range2i clip_range = Range2i::BuildWithSize(kClipRangeMinPoint, kSize);
  Array2D<Color4f> image(kSize);
  Array2D<Color4f> clipped_image = ClipArray(image, clip_range);
  EXPECT_EQ(0, clipped_image.Width());
  EXPECT_EQ(0, clipped_image.Height());
}

// Clips an array by a range that abuts the extent of the array.
TEST(Array2dUtilTest, ClipByAbuttingRange) {
  const Vector2i kSize(16, 8);
  const Point2i kClipRangeMinPoint = Point2i::Zero() + kSize;

  Range2i clip_range = Range2i::BuildWithSize(kClipRangeMinPoint, kSize);
  Array2D<Color4f> image(kSize);
  Array2D<Color4f> clipped_image = ClipArray(image, clip_range);
  EXPECT_EQ(0, clipped_image.Width());
  EXPECT_EQ(0, clipped_image.Height());
}

// Clips an array by a range of the same size as the array but offset in both x
// and y directions.
TEST(Array2dUtilTest, ClipByIntersectingRange) {
  const Vector2i kSize(16, 8);
  const Vector2i kOffByOne(1, 1);
  constexpr float kEpsilon = 1.0e-3f;
  const Point2i kClipRangeMinPoint(2, 1);

  Array2D<Color4f> image(kSize);
  for (int y = 0; y < kSize[1]; ++y) {
    for (int x = 0; x < kSize[0]; ++x) {
      float red = y * kSize[0] + x;
      image.At(x, y) = {red, 0.0f, 0.0f, 1.0f};
    }
  }

  Range2i clip_range =
      Range2i::BuildWithSize(kClipRangeMinPoint, image.GetSize() - kOffByOne);

  Array2D<Color4f> clipped_image = ClipArray(image, clip_range);
  EXPECT_EQ(14, clipped_image.Width());
  EXPECT_EQ(7, clipped_image.Height());
  for (int y = 0; y < clipped_image.Height(); ++y) {
    for (int x = 0; x < clipped_image.Width(); ++x) {
      float expected_red = (y + clip_range.GetMinPoint()[1]) * kSize[0] + x +
                           clip_range.GetMinPoint()[0];
      Color4f expected_pixel(expected_red, 0.0f, 0.0f, 1.0f);
      EXPECT_VECTOR_NEAR(expected_pixel, clipped_image.At(x, y), kEpsilon);
    }
  }
}

// Returns an image where each pixel stores its (x, y) coordinates.
Array2D<Color2i> MakePixelCoordinateImage(int row_count = 3) {
  Array2D<Color2i> image(4, row_count);
  SpatialFillArray(&image, [](const Point2i& position) {
    return Color2i(position[0], position[1]);
  });
  return image;
}

TEST(Array2dUtilTest, Transpose) {
  Array2D<Color2i> image1 = MakePixelCoordinateImage();

  Array2D<Color2i> image2;
  image2.Resize(image1.Height(), image1.Width());
  TransposeArray(image1, &image2);

  EXPECT_EQ(image1.GetSize()[0], image2.GetSize()[1]);
  EXPECT_EQ(image1.GetSize()[1], image2.GetSize()[0]);

  for (int y = 0; y < image2.Height(); ++y) {
    for (int x = 0; x < image2.Width(); ++x) {
      EXPECT_EQ(Color2i(y, x), image2.At(x, y));
    }
  }
}

TEST(Array2dUtilTest, FlipHorizontal) {
  Array2D<Color2i> image1 = MakePixelCoordinateImage();

  Array2D<Color2i> image2;
  image2.Resize(image1.GetSize());
  FlipArrayHorizontal(image1, &image2);

  EXPECT_EQ(image1.GetSize(), image2.GetSize());

  for (int y = 0; y < image2.Height(); ++y) {
    for (int x = 0; x < image2.Width(); ++x) {
      EXPECT_EQ(Color2i(image1.Width() - 1 - x, y), image2.At(x, y));
    }
  }
}

TEST(Array2dUtilTest, FlipVertical) {
  const std::array<int, 2> even_and_odd_row_counts{{3, 4}};
  for (int row_count : even_and_odd_row_counts) {
    Array2D<Color2i> image1 = MakePixelCoordinateImage(row_count);

    Array2D<Color2i> image2;
    image2.Resize(image1.GetSize());
    FlipArrayVertical(image1, &image2);

    EXPECT_EQ(image1.GetSize(), image2.GetSize());

    for (int y = 0; y < image2.Height(); ++y) {
      for (int x = 0; x < image2.Width(); ++x) {
        EXPECT_EQ(Color2i(x, image1.Height() - 1 - y), image2.At(x, y))
            << Point2i(x, y);
      }
    }

    // For odd row count, verify the middle row is unchanged.
    if (row_count % 2) {
      for (int x = 0, middle_row = row_count / 2; x < image2.Width(); ++x) {
        EXPECT_EQ(image1.At(x, middle_row), image2.At(x, middle_row));
      }
    }

    // Flip |image2| in place, and it should be equal to the original.
    FlipArrayVertical(&image2);
    EXPECT_EQ(image1.GetSize(), image2.GetSize());
    EXPECT_EQ(image1, image2);
  }
}

// Exercise SpatialForEachEntry and SpatialFill.
TEST(Array2dUtilTest, ExerciseSpatialOperations) {
  Array2D<int> input(3, 2, 0);

  // Verify initially all zero.
  SpatialForEachArrayEntry(input, [](const Point2i& position, int entry_value) {
    EXPECT_EQ(entry_value, 0) << position;
  });

  // Note that this variable is captured mutably by, and actually modified
  // inside, the lambda. This makes some assumptions about the operation of the
  // SpatialFillArray, e.g. no multi threading.
  int operation_invocations = 0;
  // Number the entries with the sum of the coordinates.
  SpatialFillArray(&input, [&operation_invocations](const Point2i& position) {
    ++operation_invocations;
    return position[0] + position[1];
  });
  EXPECT_EQ(operation_invocations, input.GetSize()[0] * input.GetSize()[1]);

  // Failure here (e.g. with zero image dimensions) would mean we aren't
  // actually exercising the EXPECT_EQ statements, in which the test could
  // otherwise return a vacuous positive test result.
  EXPECT_EQ(input.GetSize(), Vector2i(3, 2));
  // Verify the operation filled the array.
  for (int y = 0; y < input.Height(); ++y) {
    for (int x = 0; x < input.Width(); ++x) {
      EXPECT_EQ(input.At(x, y), x + y) << Point2i(x, y);
    }
  }
}

// Exercise TransformArray.
TEST(Array2dUtilTest, ExerciseTransformArray) {
  Array2D<int> input(2, 3, 0);

  // Number the entries.
  const int input_width = input.Width();
  SpatialFillArray(&input, [input_width](const Point2i& position) {
    return position[0] + position[1] * input_width;
  });
  TransformArray(input, &input,
                 [](int entry_value) { return entry_value * entry_value; });

  // Verify the operation filled the array.
  SpatialForEachArrayEntry(
      input, [&input](const Point2i& position, int entry_value) {
        const int location_number = position[0] + position[1] * input.Width();
        EXPECT_EQ(entry_value, location_number * location_number) << position;
      });
}

// Exercise SpatialTransformArray.
TEST(Array2dUtilTest, ExerciseSpatialTransformArray) {
  Array2D<int> input(2, 3, 0);

  // Number the entries.
  const int input_width = input.Width();
  SpatialFillArray(&input, [input_width](const Point2i& position) {
    return position[0] + position[1] * input_width;
  });
  SpatialTransformArray(input, &input,
                        [](const Point2i& position, int entry_value) {
                          return entry_value * entry_value * (position[0] + 1);
                        });

  // Verify the operation touched the array.
  for (int y = 0; y < input.Height(); ++y) {
    for (int x = 0; x < input.Width(); ++x) {
      const int location_number = x + y * input.Width();
      EXPECT_EQ(input.At(x, y), location_number * location_number * (x + 1))
          << Point2i(x, y);
    }
  }
}

TEST(Array2dUtilTest, TransformImage) {
  Array2D<Color2i> image1 = MakePixelCoordinateImage();

  Array2D<float> image2;

  image2.Resize(image1.GetSize());
  TransformArray(image1, &image2, [](const Color2i& p) {
    return static_cast<float>(p[0] + p[1]);
  });

  EXPECT_EQ(image1.GetSize(), image2.GetSize());

  for (int y = 0; y < image2.Height(); ++y) {
    for (int x = 0; x < image2.Width(); ++x) {
      EXPECT_EQ(static_cast<float>(image1.At(x, y)[0] + image1.At(x, y)[1]),
                image2.At(x, y));
    }
  }
}

TEST(Array2dUtilTest, IntegralImage_Uniform) {
  // Test with a uniform image.
  Array2D<Color2i> image(6, 5);
  image.Fill({3, -2});

  Array2D<Color2i> integral;

  integral.Resize(image.GetSize());
  ComputeIntegralArray(image, &integral);

  EXPECT_EQ(image.GetSize(), integral.GetSize());

  for (int y = 0; y < image.Height(); ++y) {
    for (int x = 0; x < image.Width(); ++x) {
      EXPECT_EQ((x + 1) * (y + 1) * Color2i(3, -2), integral.At(x, y));
    }
  }

  EXPECT_EQ(Color2i(3, -2), SumFromIntegralArray(integral, {{0, 0}, {0, 0}}));

  EXPECT_EQ(Color2i(3, -2) * 4,
            SumFromIntegralArray(integral, {{3, 2}, {4, 3}}));
}

TEST(Array2dUtilTest, IntegralImage_Specific) {
  // Test with a specific example image.
  std::vector<int> image_data = {
      1,  0,  -3, 1,  //
      -2, -3, 0,  1,  //
      3,  -1, -2, 3,  //
      -1, 1,  2,  4,  //
  };

  Array2D<int> image(4, 4);
  std::copy(image_data.begin(), image_data.end(), image.Data());

  Array2D<int> integral;
  integral.Resize(image.GetSize());
  ComputeIntegralArray(image, &integral);

  // 2x2 square in the middle.
  EXPECT_EQ(-3 + 0 - 1 - 2, SumFromIntegralArray(integral, {{1, 1}, {2, 2}}));

  // Top edge.
  EXPECT_EQ(1 + 0 - 3 + 1, SumFromIntegralArray(integral, {{0, 0}, {3, 0}}));

  // Right edge.
  EXPECT_EQ(1 + 1 + 3 + 4, SumFromIntegralArray(integral, {{3, 0}, {3, 3}}));

  // Left edge.
  EXPECT_EQ(1 - 2 + 3 - 1, SumFromIntegralArray(integral, {{0, 0}, {0, 3}}));

  // Bottom edge.
  EXPECT_EQ(-1 + 1 + 2 + 4, SumFromIntegralArray(integral, {{0, 3}, {3, 3}}));

  // Top-left 3x3 subregion.
  EXPECT_EQ(0 +                 //
                (1 + 0 - 3) +   //
                (-2 - 3 + 0) +  //
                (3 - 1 - 2),    //
            SumFromIntegralArray(integral, {{0, 0}, {2, 2}}));
}

TEST(Array2dUtil, BorderExtendedAt) {
  constexpr float kEpsilon = 1.0e-5f;
  // clang-format off
  std::vector<float> array_data = { 1.0f,  2.0f,  3.0f,  4.0f,
                                    5.0f,  6.0f,  7.0f,  8.0f,
                                    9.0f, 10.0f, 11.0f, 12.0f };
  // clang-format on
  Array2D<float> array(4, 3);
  std::copy(array_data.begin(), array_data.end(), array.Data());

  // clang-format off
  EXPECT_NEAR(10.0f, BorderExtendedAt(array,  1,  2), kEpsilon);
  EXPECT_NEAR( 9.0f, BorderExtendedAt(array, -1,  2), kEpsilon);
  EXPECT_NEAR( 1.0f, BorderExtendedAt(array, -5, -3), kEpsilon);
  EXPECT_NEAR( 2.0f, BorderExtendedAt(array,  1, -1), kEpsilon);
  EXPECT_NEAR( 9.0f, BorderExtendedAt(array, -1,  6), kEpsilon);
  EXPECT_NEAR( 9.0f, BorderExtendedAt(array, -2,  3), kEpsilon);
  EXPECT_NEAR( 1.0f, BorderExtendedAt(array, -1, -2), kEpsilon);
  EXPECT_NEAR(12.0f, BorderExtendedAt(array,  3,  4), kEpsilon);
  EXPECT_NEAR( 9.0f, BorderExtendedAt(array, -3,  5), kEpsilon);
  EXPECT_NEAR(12.0f, BorderExtendedAt(array,  4,  5), kEpsilon);
  // clang-format on
}

TEST(Array2dUtil, ComputeRowAverage) {
  constexpr float kEpsilon = 1.0e-5f;
  Array2D<float> array(6, 5, 2.0f);

  EXPECT_NEAR(2.0f, ComputeRowAverage(array, Range1f(0.3, 0.9), 2), kEpsilon);

  EXPECT_NEAR(2.0f, ComputeRowAverage(array, Range1f(0.8, 4.1), 3), kEpsilon);
}

TEST(Array2dUtil, ComputeBoxAverage) {
  constexpr float kEpsilon = 1.0e-5f;
  std::vector<float> array_data = {2.0f, 3.0f, 0.0f, 1.0f,  //
                                   4.0f, 1.0f, 2.0f, 0.0f,  //
                                   3.0f, 5.0f, 1.0f, 2.0f};
  Array2D<float> array(4, 3);
  std::copy(array_data.begin(), array_data.end(), array.Data());

  EXPECT_NEAR(2.0f, ComputeBoxAverage(array, {-0.5f, 3.5f}, {-0.5f, 2.5f}),
              kEpsilon);

  EXPECT_NEAR(3.4375f, ComputeBoxAverage(array, {-0.3f, 0.4f}, {0.2f, 1.8f}),
              kEpsilon);
}

}  // namespace
}  // namespace base
}  // namespace seurat
