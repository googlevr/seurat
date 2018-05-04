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

// Array Utilities are STL-style generic algorithms on 2D array containers,
// particularly base::Array2D.
//
// All versions of algorithms that are spatially varying include the word
// spatial in the name. The coordinate-independent versions have simpler names.

#ifndef VR_SEURAT_BASE_ARRAY2D_UTIL_H_
#define VR_SEURAT_BASE_ARRAY2D_UTIL_H_

#include <utility>

#include "ion/math/range.h"
#include "ion/math/rangeutils.h"
#include "ion/math/vector.h"
#include "seurat/base/array2d.h"

namespace seurat {
namespace base {

// Fills entries in |result| with the value of a given |function| of the entry's
// coordinates.
//
// result_x,y = f(x, y) for all (x,y) in result's domain.
//
// Implements an analog of std::fill over the array.
template <typename ArrayType, typename UnarySpatialOp>
void SpatialFillArray(ArrayType* result, const UnarySpatialOp function) {
  for (int y = 0; y < result->Height(); ++y) {
    for (int x = 0; x < result->Width(); ++x) {
      result->At(x, y) = function(ion::math::Point2i(x, y));
    }
  }
}

// Evaluates a coordinate-dependent |function| on all of the entries in |array|.
//
// Evaluates f((x, y), array_x,y) for all (x,y) in array's domain.
//
// Implements an analog of std::fill over the array.
template <typename ArrayType, typename UnarySpatialOp>
void SpatialForEachArrayEntry(const ArrayType& array,
                              const UnarySpatialOp function) {
  for (int y = 0; y < array.Height(); ++y) {
    for (int x = 0; x < array.Width(); ++x) {
      function(ion::math::Point2i(x, y), array.At(x, y));
    }
  }
}

// Copies elements from the |src| array into the |dst| array at the specified
// |dst_offset|.
//
// Note:  No clipping is performed.  The specified |src| and offset must fit
// within the |dst| array.
template <typename SrcArrayType, typename DstArrayType>
void CopyArray(const SrcArrayType& src, DstArrayType* dst,
               const ion::math::Vector2i& dst_offset) {
  CHECK_LE(src.Width() + dst_offset[0], dst->Width());
  CHECK_LE(src.Height() + dst_offset[1], dst->Height());
  CHECK_GE(dst_offset[0], 0);
  CHECK_GE(dst_offset[1], 0);
  for (int y = 0; y < src.Height(); ++y) {
    for (int x = 0; x < src.Width(); ++x) {
      dst->At(x + dst_offset[0], y + dst_offset[1]) = src.At(x, y);
    }
  }
}

// Returns the given |source_array| clipped to the |clip_range|. Part of
// |clip_range| may be outside the extent of |source_array|.
template <typename ArrayType>
ArrayType ClipArray(const ArrayType& source_array,
                    const ion::math::Range2i& clip_range) {
  // The size of an ion::math::Range2i is smaller by one in each dimension than
  // the size of the corresponding Array2D.
  const ion::math::Vector2i kOffByOne(1, 1);
  ion::math::Range2i array_range = ion::math::Range2i::BuildWithSize(
      ion::math::Point2i::Zero(), source_array.GetSize() - kOffByOne);
  ion::math::Range2i result_range =
      ion::math::RangeIntersection(array_range, clip_range);

  // Returns an empty array if the resulting range is empty.
  if (result_range.IsEmpty()) {
    return ArrayType();
  }

  // Returns the |source_array| when there is no clipping.
  if (result_range.GetSize() + kOffByOne == source_array.GetSize()) {
    return source_array;
  }

  ArrayType result(result_range.GetSize() + kOffByOne);
  ion::math::Point2i result_origin = result_range.GetMinPoint();
  for (int y = 0; y < result.Height(); ++y) {
    for (int x = 0; x < result.Width(); ++x) {
      ion::math::Point2i element_position(x, y);
      result.At(element_position) =
          source_array.At(result_origin + element_position);
    }
  }
  return result;
}

// Swaps the (x, y) elements in |original| into (y,x ) elements in |transposed|.
// Automatically resizes the destination |transposed| to match the (height,
// width) of |original|.
template <typename ArrayType>
void TransposeArray(const ArrayType& original, ArrayType* transposed) {
  transposed->Resize(original.Height(), original.Width());
  for (int y = 0; y < original.Height(); ++y) {
    for (int x = 0; x < original.Width(); ++x) {
      transposed->At(y, x) = original.At(x, y);
    }
  }
}

// Flips |original| horizontally into |flipped|. Automatically resizes the
// destination |flipped| to match the size of the |original|.
template <typename ArrayType>
void FlipArrayHorizontal(const ArrayType& original, ArrayType* flipped) {
  flipped->Resize(original.Width(), original.Height());
  for (int y = 0; y < original.Height(); ++y) {
    for (int x = 0; x < original.Width(); ++x) {
      flipped->At(x, y) = original.At(original.Width() - x - 1, y);
    }
  }
}

// Flips |flipped_in_place| vertically into itself.
template <typename ArrayType>
void FlipArrayVertical(ArrayType* flipped_in_place) {
  for (int y = 0; y < flipped_in_place->Height() / 2; ++y) {
    for (int x = 0; x < flipped_in_place->Width(); ++x) {
      auto& low_element = flipped_in_place->At(x, y);
      auto& high_element =
          flipped_in_place->At(x, flipped_in_place->Height() - 1 - y);
      using std::swap;
      swap(low_element, high_element);
    }
  }
}

// Flips |original| vertically into |flipped|. Automatically resizes the
// destination |flipped| to match the size of the |original|.
template <typename ArrayType>
void FlipArrayVertical(const ArrayType& original, ArrayType* flipped) {
  flipped->Resize(original.Width(), original.Height());
  for (int y = 0; y < original.Height(); ++y) {
    for (int x = 0; x < original.Width(); ++x) {
      flipped->At(x, y) = original.At(x, original.Height() - y - 1);
    }
  }
}

// Evaluates a given |function| (independent of coordinates) on each element of
// the |original| array and stores the output in |result|. Automatically resizes
// the destination |result| to match the size of the |original|.
//
// result_x,y = f(original_x,y) for all (x,y) in result's domain.
//
// Implements an analog of std::transform over the array.
template <typename ArrayTypeIn, typename ArrayTypeOut, typename UnaryOp>
void TransformArray(const ArrayTypeIn& original, ArrayTypeOut* result,
                    const UnaryOp function) {
  result->Resize(original.Width(), original.Height());
  for (int y = 0; y < original.Height(); ++y) {
    for (int x = 0; x < original.Width(); ++x) {
      result->At(x, y) = function(original.At(x, y));
    }
  }
}

// Evaluates a given spatial unary |operation| on each element of the |original|
// array and stores the output in |result|. Automatically resizes the
// destination |result| to match the size of the |original|.
//
// result_x,y = f((x, y), original_x,y) for all (x,y) in result's  domain.
//
// Implements an analog of std::transform over the array.
template <typename ArrayTypeIn, typename ArrayTypeOut, typename UnarySpatialOp>
void SpatialTransformArray(const ArrayTypeIn& original, ArrayTypeOut* result,
                           const UnarySpatialOp operation) {
  result->Resize(original.Width(), original.Height());
  for (int y = 0; y < original.Height(); ++y) {
    for (int x = 0; x < original.Width(); ++x) {
      result->At(x, y) = operation(ion::math::Point2i(x, y), original.At(x, y));
    }
  }
}

// Computes an integral array, where each element, at (x, y), is the sum of all
// elements of the |original| array which are in the subregion from (0,0) to
// (x, y), inclusive.
template <typename ArrayType>
void ComputeIntegralArray(const ArrayType& original,
                          ArrayType* integral_array) {
  integral_array->Resize(original.Width(), original.Height());
  for (int y = 0; y < original.Height(); ++y) {
    for (int x = 0; x < original.Width(); ++x) {
      auto sum = original.At(x, y);
      if (y > 0) {
        sum += integral_array->At(x, y - 1);
      }
      if (x > 0) {
        sum += integral_array->At(x - 1, y);
      }
      if (x > 0 && y > 0) {
        sum -= integral_array->At(x - 1, y - 1);
      }
      integral_array->At(x, y) = sum;
    }
  }
}

// Given an |integral_array| computes the sum of the original array over the
// specified |element_range| which is *inclusive* on both ends.
template <typename ArrayType>
typename ArrayType::ElementType SumFromIntegralArray(
    const ArrayType& integral_array, const ion::math::Range2i element_range) {
  ion::math::Point2i corner = element_range.GetMaxPoint();
  typename ArrayType::ElementType partial_sum = integral_array.At(corner);

  corner = {element_range.GetMaxPoint()[0], element_range.GetMinPoint()[1] - 1};
  if (integral_array.IsInside(corner)) {
    partial_sum -= integral_array.At(corner);
  }

  corner = {element_range.GetMinPoint()[0] - 1, element_range.GetMaxPoint()[1]};
  if (integral_array.IsInside(corner)) {
    partial_sum -= integral_array.At(corner);
  }

  corner = {element_range.GetMinPoint()[0] - 1,
            element_range.GetMinPoint()[1] - 1};
  if (integral_array.IsInside(corner)) {
    partial_sum += integral_array.At(corner);
  }

  return partial_sum;
}

// Returns the |array| element at indices x and y, when both indices are within
// their valid ranges. If one of, or both, x and y point outside the array, this
// function returns the element on the border of the array that is closest to
// the location (x,y) requested. This is equivalent to extending the array by
// replicating border elements.
template <typename ArrayType>
typename ArrayType::ElementType BorderExtendedAt(const ArrayType& array, int x,
                                                 int y) {
  if (x < 0) {
    return BorderExtendedAt(array, 0, y);
  }
  if (array.Width() <= x) {
    return BorderExtendedAt(array, array.Width() - 1, y);
  }
  if (y < 0) {
    return BorderExtendedAt(array, x, 0);
  }
  if (array.Height() <= y) {
    return BorderExtendedAt(array, x, array.Height() - 1);
  }
  return array.At(x, y);
}

// Returns the average of array elements within the |y| row and within range
// |x_range|. The endpoints of the x range need not be integers: partially
// covered entries will be weighted linearly. The samples being averaged are
// assumed to be centered at integer x and y.
template <typename ArrayType>
typename ArrayType::ElementType ComputeRowAverage(
    const ArrayType& array, const ion::math::Range1f& x_range, int y) {
  // Lowest-x sample that can contribute to the average over x_range.
  const int min_x = static_cast<int>(std::floor(x_range.GetMinPoint() + 0.5f));
  // Highest-x sample that can contribute to the average over x_range.
  const int max_x = static_cast<int>(std::floor(x_range.GetMaxPoint() + 0.5f));

  if (min_x == max_x) {
    // The x range spans one single entry.
    return BorderExtendedAt(array, min_x, y);
  }

  typename ArrayType::ElementType average;

  // Accumulate the first partially covered entry.
  const float weight_first_partial_pixel = min_x + 0.5f - x_range.GetMinPoint();
  average = weight_first_partial_pixel * BorderExtendedAt(array, min_x, y);

  // Accumulate fully covered entries.
  for (int x = 1 + min_x; x < max_x; ++x) {
    average += BorderExtendedAt(array, x, y);
  }

  // Accumulate the last partially covered entry.
  const float weight_last_partial_pixel = x_range.GetMaxPoint() + 0.5f - max_x;
  average += weight_last_partial_pixel * BorderExtendedAt(array, max_x, y);

  // Normalize
  average /= x_range.GetSize();

  return average;
}

// Returns the average of array elements within the ranges |x_range|
// and |y_range|. The endpoints of the x and y ranges need not be integers:
// partially covered entries will be weighted bilinearly.
template <typename ArrayType>
typename ArrayType::ElementType ComputeBoxAverage(
    const ArrayType& array, const ion::math::Range1f& x_range,
    const ion::math::Range1f& y_range) {
  const int min_y = static_cast<int>(std::floor(y_range.GetMinPoint() + 0.5f));
  const int max_y = static_cast<int>(std::floor(y_range.GetMaxPoint() + 0.5f));

  if (min_y == max_y) {
    // The y range spans one single row.
    return ComputeRowAverage(array, x_range, min_y);
  }

  typename ArrayType::ElementType average;

  // Accumulate first partially covered rows.
  const float weight_first_partial_line = min_y + 0.5f - y_range.GetMinPoint();
  average =
      weight_first_partial_line * ComputeRowAverage(array, x_range, min_y);

  // Accumulate fully covered rows.
  for (int i = 1 + min_y; i < max_y; ++i) {
    average += ComputeRowAverage(array, x_range, i);
  }

  // Accumulate last partially covered row.
  const float weight_last_partial_line = y_range.GetMaxPoint() + 0.5f - max_y;
  average +=
      weight_last_partial_line * ComputeRowAverage(array, x_range, max_y);

  // Normalize.
  average /= y_range.GetSize();

  return average;
}

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_ARRAY2D_UTIL_H_
