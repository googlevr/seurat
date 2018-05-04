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

#include "seurat/image/filter.h"
#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "gtest/gtest.h"

using ion::math::Point2f;
using ion::math::Range2f;
using ion::math::Vector2f;

namespace seurat {
namespace image {

TEST(Filter, BoxFilter) {
  BoxFilter filter(Vector2f(1.0f, 1.0f));
  const float expected_radius = std::sqrt(0.5f * 0.5f + 0.5f * 0.5f);
  const Range2f expected_range =
      Range2f(Point2f(-0.5f, -0.5f), Point2f(0.5f, 0.5f));
  EXPECT_EQ(expected_radius, filter.GetRadius());
  EXPECT_EQ(expected_range, filter.GetRange());
  EXPECT_EQ(1.0f, filter.Eval(Point2f(0.0f, 0.0f)));
  EXPECT_EQ(1.0f, filter.Eval(Point2f(-0.5f, -0.5f)));
  EXPECT_EQ(1.0f, filter.Eval(Point2f(0.5f, -0.5f)));
  EXPECT_EQ(1.0f, filter.Eval(Point2f(0.5f, 0.5f)));
  EXPECT_EQ(1.0f, filter.Eval(Point2f(-0.5f, 0.5f)));
  EXPECT_EQ(0.0f, filter.Eval(Point2f(0.50001f, 0.0f)));
  EXPECT_EQ(0.0f, filter.Eval(Point2f(-0.50001f, 0.0f)));
  EXPECT_EQ(0.0f, filter.Eval(Point2f(0.0f, 0.50001f)));
  EXPECT_EQ(0.0f, filter.Eval(Point2f(0.0f, -0.50001f)));
}

TEST(Filter, BSplineFilter) {
  BSplineFilter filter;
  const float expected_radius = 2.0f;
  const Range2f expected_range =
      Range2f(Point2f(-2.0f, -2.0f), Point2f(2.0f, 2.0f));
  EXPECT_EQ(expected_radius, filter.GetRadius());
  EXPECT_EQ(expected_range, filter.GetRange());
  EXPECT_EQ(0.0f, filter.Eval(Point2f(0.0f, 2.0f)));
  EXPECT_EQ(0.0f, filter.Eval(Point2f(0.0f, -2.0f)));
  EXPECT_EQ(0.0f, filter.Eval(Point2f(2.0f, 0.0f)));
  EXPECT_EQ(0.0f, filter.Eval(Point2f(-2.0f, 0.0f)));
  EXPECT_EQ(4.0f / 6.0f, filter.Eval(Point2f(0.0f, 0.0f)));
  EXPECT_EQ(1.0f / 6.0f, filter.Eval(Point2f(0.0f, 1.0f)));
  EXPECT_EQ(1.0f / 6.0f, filter.Eval(Point2f(-1.0f, 0.0f)));
}

TEST(Filter, MitchellFilter) {
  MitchellFilter filter;
  const float expected_radius = 2.0f;
  const Range2f expected_range =
      Range2f(Point2f(-2.0f, -2.0f), Point2f(2.0f, 2.0f));
  EXPECT_EQ(expected_radius, filter.GetRadius());
  EXPECT_EQ(expected_range, filter.GetRange());
  EXPECT_EQ(0.0f, filter.Eval(Point2f(0.0f, 2.0f)));
  EXPECT_EQ(0.0f, filter.Eval(Point2f(0.0f, -2.0f)));
  EXPECT_EQ(0.0f, filter.Eval(Point2f(2.0f, 0.0f)));
  EXPECT_EQ(0.0f, filter.Eval(Point2f(-2.0f, 0.0f)));
  EXPECT_EQ((6.0f - 2.0f / 3.0f) / 6.0f, filter.Eval(Point2f(0.0f, 0.0f)));
  EXPECT_NEAR(1.0f / 18.0f, filter.Eval(Point2f(0.0f, 1.0f)), 1e-5f);
  EXPECT_NEAR(1.0f / 18.0f, filter.Eval(Point2f(-1.0f, 0.0f)), 1e-5f);
}

TEST(Filter, GaussianFilter) {
  const float kSigma = 1.0f;
  const float kRadius = 3.0f;
  GaussianFilter filter(kSigma, kRadius);
  EXPECT_EQ(kRadius, filter.GetRadius());
  const Range2f expected_range(Point2f(-kRadius, -kRadius),
                               Point2f(kRadius, kRadius));
  EXPECT_EQ(expected_range, filter.GetRange());
  EXPECT_EQ(0.0f, filter.Eval(Point2f(0.0f, 3.0f)));
  EXPECT_EQ(0.0f, filter.Eval(Point2f(0.0f, 4.0f)));
  EXPECT_EQ(0.0f, filter.Eval(Point2f(-3.0f, 0.0f)));
  EXPECT_EQ(0.0f, filter.Eval(Point2f(-4.0f, 0.0f)));
  const float value_at_cutoff =
      1.0f / std::sqrt(2.0f * M_PI) *
      std::exp(-kRadius * kRadius / (2.0 * kSigma * kSigma));
  const float value_at_center = 1.0f / std::sqrt(2.0f * M_PI) - value_at_cutoff;
  EXPECT_FLOAT_EQ(value_at_center, filter.Eval(Point2f(0.0f, 0.0f)));
}

TEST(Filter, Wendland31Filter) {
  const float kRadius = 2.5f;
  Wendland31Filter filter(kRadius);
  EXPECT_EQ(kRadius, filter.GetRadius());
  const Range2f expected_range(Point2f(-kRadius, -kRadius),
                               Point2f(kRadius, kRadius));
  EXPECT_EQ(expected_range, filter.GetRange());
  std::vector<Point2f> outside_support_points = {
      {0.0f, 3.0f}, {0.0f, 2.5f}, {-2.5f, 0.0f}, {-3.0f, 0.0f}};
  for (const auto& p : outside_support_points) EXPECT_EQ(0.0f, filter.Eval(p));
  EXPECT_EQ(1.0f, filter.Eval(Point2f(0.0f, 0.0f)));
  float r = 1.0f / kRadius;
  float omr = 1.0f - r;
  float omr2 = omr * omr;
  float val = omr2 * omr2 * (1.0f + 4.0f * r);
  EXPECT_NEAR(val, filter.Eval(Point2f(0.0f, 1.0f)), 1e-5f);
  EXPECT_NEAR(val, filter.Eval(Point2f(-1.0f, 0.0f)), 1e-5f);
  EXPECT_NEAR(val, filter.Eval(Point2f(0.5f, 0.5f * std::sqrt(3.0f))), 1e-5f);
}

}  // namespace image
}  // namespace seurat
