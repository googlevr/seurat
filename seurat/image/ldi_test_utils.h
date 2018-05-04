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

#ifndef VR_SEURAT_IMAGE_LDI_TEST_UTILS_H_
#define VR_SEURAT_IMAGE_LDI_TEST_UTILS_H_

#include <algorithm>
#include <functional>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/image/ldi.h"

namespace seurat {
namespace image {

using ::testing::AssertionFailure;
using ::testing::AssertionResult;
using ::testing::AssertionSuccess;
using ::testing::PrintToString;
using ion::math::Vector2i;

template <typename T>
AssertionResult LdiValid(const Ldi<T>& ldi) {
  const Vector2i size = ldi.GetSize();
  const int sample_count = ldi.GetSampleCount();
  int64 actual_sample_count = 0;
  for (int y = 0; y < size[1]; ++y) {
    for (int x = 0; x < size[0]; ++x) {
      const int pixel_sample_count = ldi.GetSampleCount({x, y});
      if (pixel_sample_count < 0) {
        return AssertionFailure()
               << "GetSampleCount({" << x << ", " << y
               << "}) must be non-negative: value is " << pixel_sample_count;
      }
      actual_sample_count += pixel_sample_count;
    }
  }
  if (actual_sample_count != sample_count) {
    return AssertionFailure()
           << "GetSampleCount() does not match actual_sample_count: "
           << sample_count << " vs. " << actual_sample_count;
  }

  return AssertionSuccess();
}

template <typename T, typename ColorEqual = std::equal_to<T>,
          typename DepthEqual = std::equal_to<typename T::ValueType>>
AssertionResult LdiEquals(const Ldi<T>& lhs, const Ldi<T>& rhs,
                          ColorEqual color_equal = ColorEqual(),
                          DepthEqual depth_equal = DepthEqual()) {
  AssertionResult result = LdiValid(lhs);
  if (!result) {
    return AssertionFailure() << "lhs is invalid: " << result.message();
  }
  result = LdiValid(rhs);
  if (!result) {
    return AssertionFailure() << "rhs is invalid: " << result.message();
  }

  if (lhs.GetSize() != rhs.GetSize()) {
    return AssertionFailure()
           << "lhs.GetSize() does not match rhs.GetSize(): " << lhs.GetSize()
           << " vs. " << rhs.GetSize();
  }
  if (lhs.GetSampleCount() != rhs.GetSampleCount()) {
    return AssertionFailure()
           << "lhs.GetSampleCount() does not match rhs.GetSampleCount(): "
           << lhs.GetSampleCount() << " vs. " << rhs.GetSampleCount();
  }

  for (int y = 0; y < lhs.GetHeight(); ++y) {
    for (int x = 0; x < lhs.GetWidth(); ++x) {
      const auto& lhs_colors = lhs.GetColors({x, y});
      const auto& rhs_colors = rhs.GetColors({x, y});
      if (lhs_colors.size() != rhs_colors.size() ||
          !std::equal(lhs_colors.begin(), lhs_colors.end(), rhs_colors.begin(),
                      color_equal)) {
        return AssertionFailure() << "lhs.GetColors({" << x << ", " << y
                                  << "}) does not match rhs.GetColors({" << x
                                  << ", " << y << "}):\n"
                                  << PrintToString(lhs_colors) << "\n  vs.\n"
                                  << PrintToString(rhs_colors);
      }
      const auto& lhs_depths = lhs.GetDepths({x, y});
      const auto& rhs_depths = rhs.GetDepths({x, y});
      if (lhs_depths.size() != rhs_depths.size() ||
          !std::equal(lhs_depths.begin(), lhs_depths.end(), rhs_depths.begin(),
                      depth_equal)) {
        return AssertionFailure() << "lhs.GetDepths({" << x << ", " << y
                                  << "}) does not match rhs.GetDepths({" << x
                                  << ", " << y << "}):\n"
                                  << PrintToString(lhs_depths) << "\n  vs.\n"
                                  << PrintToString(rhs_depths);
      }
    }
  }

  return AssertionSuccess();
}

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_LDI_TEST_UTILS_H_
