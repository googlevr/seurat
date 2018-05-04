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

#ifndef VR_SEURAT_TESTING_ION_TEST_UTILS_H_
#define VR_SEURAT_TESTING_ION_TEST_UTILS_H_

#include <numeric>

#include "ion/math/vector.h"
#include "gtest/gtest.h"

namespace seurat {
namespace testing {

// VectorNear is a 3-argument predicate formatter function overload for
// EXPECT_PRED_FORMAT3, which supports analogous behavior to EXPECT_NEAR for
// vector types.
//
// See "Using a Predicate-Formatter" in Advanced Topics:
// https://github.com/google/googletest/blob/master/googletest/docs/AdvancedGuide.md#using-a-predicate-formatter
template <typename T>
::testing::AssertionResult VectorNear(const std::string& expr1,
                                      const std::string& expr2,
                                      const std::string& abs_error_expr,
                                      const T& val1, const T& val2,
                                      typename T::ValueType abs_error) {
  for (int i = 0; i < T::kDimension; ++i) {
    if (std::abs(val1[i] - val2[i]) > abs_error) {
      return ::testing::AssertionFailure()
             << "The difference between " << expr1 << " and " << expr2 << " is "
             << (val1 - val2) << ", which exceeds " << abs_error_expr << " at ["
             << i << "], where\n"
             << expr1 << " evaluates to " << val1 << ",\n"
             << expr2 << " evaluates to " << val2 << ", and\n"
             << abs_error_expr << " evaluates to " << abs_error << ".";
    }
  }

  return ::testing::AssertionSuccess();
}

// Similar to EXPECT_NEAR, but for Ion Vector and Point types.
#define EXPECT_VECTOR_NEAR(val1, val2, abs_error) \
  EXPECT_PRED_FORMAT3(::seurat::testing::VectorNear, val1, val2, abs_error)

// VectorEqual is a 3-argument predicate formatter function overload for
// EXPECT_PRED_FORMAT3, which supports analogous behavior to EXPECT_FLOAT_EQ
// for vector types. Therefore, despite including Equal in the name, it tests
// for a near-equality with a relative tolerance specified in units-last-place
// or |ulps_tolerance|.
template <typename T>
::testing::AssertionResult VectorEqual(const std::string& expr1,
                                       const std::string& expr2,
                                       const std::string& abs_error_expr,
                                       const T& val1, const T& val2,
                                       typename T::ValueType ulps_tolerance) {
  for (int i = 0; i < T::kDimension; ++i) {
    // Construct an absolute tolerance from the specified ULP error factor,
    // using the magnitude of the second (by convention, the expected) value.
    // Note: If FloatingPoint<T> were not internal testing, we'd probably use
    // that.
    const typename T::ValueType kAbsoluteError =
        std::numeric_limits<typename T::ValueType>::epsilon() * ulps_tolerance *
        std::abs(val2[i]);
    if (std::abs(val1[i] - val2[i]) > kAbsoluteError) {
      return ::testing::AssertionFailure()
             << "The difference between " << expr1 << " and " << expr2 << " is "
             << (val1 - val2) << ", which exceeds " << kAbsoluteError << " at ["
             << i << "], where\n"
             << expr1 << " evaluates to " << val1 << ",\n"
             << expr2 << " evaluates to " << val2 << ", and\n"
             << abs_error_expr << " ULPs evaluates to " << kAbsoluteError
             << " in absolute error for the first failing component.";
    }
  }

  return ::testing::AssertionSuccess();
}

// Similar to EXPECT_FLOAT_EQ, but for Ion Vector and Point types.
// Note: we pass in 4.0 ULP as the relative tolerance to match EXPECT_FLOAT_EQ.
#define EXPECT_VECTOR_FLOAT_EQ(val1, val2) \
  EXPECT_PRED_FORMAT3(::seurat::testing::VectorEqual, val1, val2, 4.0f)

}  // namespace testing
}  // namespace seurat

#endif  // VR_SEURAT_TESTING_ION_TEST_UTILS_H_
