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

#include "seurat/image/low_discrepancy_sampling.h"

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "gtest/gtest.h"

using ion::math::Point2f;
using ion::math::Range2f;

namespace seurat {
namespace image {
namespace {

TEST(LowDiscrepancySampling, LarcherPillichshammer) {
  const uint32 kLog2NumSamples = 8U;
  const uint32 kNumSamples = 1U << kLog2NumSamples;
  const uint32 kScramblerA = 0U;
  const uint32 kScramblerB = 579135792U;  // Arbitrary random number.
  // Generate all elementary intervals and check that each of them
  // contains exactly one sample from each point set. And check that
  // the two point sets are different.
  for (int x_exp = 0; x_exp <= kLog2NumSamples; ++x_exp) {
    int y_exp = kLog2NumSamples - x_exp;
    for (int x = 0; x < (1 << x_exp); ++x) {
      for (int y = 0; y < (1 << y_exp); ++y) {
        const float x_min = x * 1.0f / static_cast<float>(1 << x_exp);
        const float x_max = (x + 1) * 1.0f / static_cast<float>(1 << x_exp);
        const float y_min = y * 1.0f / static_cast<float>(1 << y_exp);
        const float y_max = (y + 1) * 1.0f / static_cast<float>(1 << y_exp);
        bool has_sample_a_in_interval = false;
        bool has_sample_b_in_interval = false;
        for (int i = 0; i < kNumSamples; ++i) {
          const Point2f sample_position_a =
              LarcherPillichshammer2D(i, kNumSamples, kScramblerA);
          const Point2f sample_position_b =
              LarcherPillichshammer2D(i, kNumSamples, kScramblerB);
          EXPECT_NE(sample_position_a, sample_position_b);
          if (sample_position_a[0] >= x_min && sample_position_a[0] < x_max &&
              sample_position_a[1] >= y_min && sample_position_a[1] < y_max) {
            EXPECT_FALSE(has_sample_a_in_interval);
            has_sample_a_in_interval = true;
          }
          if (sample_position_b[0] >= x_min && sample_position_b[0] < x_max &&
              sample_position_b[1] >= y_min && sample_position_b[1] < y_max) {
            EXPECT_FALSE(has_sample_b_in_interval);
            has_sample_b_in_interval = true;
          }
        }
        EXPECT_TRUE(has_sample_a_in_interval);
        EXPECT_TRUE(has_sample_b_in_interval);
      }
    }
  }
}

}  // namespace
}  // namespace image
}  // namespace seurat
