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

#include "seurat/testing/ion_test_utils.h"

#include <numeric>

#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {

// Exercise EXPECT_VECTOR_NEAR and VectorNear.
TEST(IonTestUtilsTest, VectorComparisonPredicate) {
  const ion::math::Vector3f kInitial(1.1f, 3.14f, -1.997f);
  const ion::math::Vector3f half_initial = kInitial * 0.5f;
  EXPECT_VECTOR_NEAR(kInitial, half_initial, 1.58f);

  EXPECT_THAT(
      seurat::testing::VectorNear("kInitial", "half_initial", "1.0f",
                                     kInitial, half_initial, 1.0f)
          .message(),
      ::testing::HasSubstr("The difference between kInitial and half_initial"));

  const ion::math::Point3f kInitialPoint(1.1f, 0.0f, 0.0f);
  const ion::math::Point3f half_initial_point = kInitialPoint * 0.5f;
  EXPECT_VECTOR_NEAR(kInitialPoint, half_initial_point, 0.55f);
}

// Exercise EXPECT_VECTOR_FLOAT_EQUAL and VectorEqual.
TEST(IonTestUtilsTest, FloatVectorEqual) {
  const ion::math::Vector3f kInitial(1.1f, 3023409.14f, -1.997f);
  EXPECT_VECTOR_FLOAT_EQ(kInitial, kInitial);

  // Test the range of values up to several units-last place around the precise
  // initial value.
  const float kFltEpsilon = std::numeric_limits<float>::epsilon();
  for (int num_ulps = -3; num_ulps < 4; ++num_ulps) {
    // Note that at +/- 4 ULPS (which this loop doesn't hit), round off makes
    // the modification actually go slightly outside 4 ULPs error.
    const ion::math::Vector3f modified_by_num_ulps =
        kInitial + kInitial * kFltEpsilon * num_ulps;
    EXPECT_VECTOR_FLOAT_EQ(modified_by_num_ulps, kInitial)
        << " At " << num_ulps;
  }

  const ion::math::Vector3f modify_initial_5ulps_error =
      kInitial + kInitial * kFltEpsilon * 5.0f;
  EXPECT_THAT(seurat::testing::VectorEqual(
                  "kInitial", "modify_initial_5ulps_error", "4.0f", kInitial,
                  modify_initial_5ulps_error, 4.0f)
                  .message(),
              ::testing::HasSubstr("in absolute error for the first"));
}

}  // namespace
