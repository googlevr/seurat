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

#include "seurat/geometry/convex_hull2d.h"

#include <vector>

#include "ion/math/vector.h"
#include "gtest/gtest.h"

using ion::math::Point2f;

namespace seurat {
namespace geometry {

TEST(ConvexHull2dTest, TestConvexHull) {
  // A set of points consisting of a unit-square and several interior points.
  std::vector<Point2f> points;
  points.push_back(Point2f(0.0f, 0.0f));
  points.push_back(Point2f(1.0f, 0.0f));
  points.push_back(Point2f(0.0f, 1.0f));
  points.push_back(Point2f(1.0f, 1.0f));

  points.push_back(Point2f(0.1f, 0.3f));
  points.push_back(Point2f(0.2f, 0.6f));
  points.push_back(Point2f(0.8f, 0.1f));
  points.push_back(Point2f(0.9f, 0.9f));
  points.push_back(Point2f(0.9f, 0.7f));

  auto hull = ComputeConvexHull(points);

  EXPECT_EQ(4, hull.size());
  EXPECT_EQ(1, std::count(hull.begin(), hull.end(), Point2f(0.0f, 0.0f)));
  EXPECT_EQ(1, std::count(hull.begin(), hull.end(), Point2f(1.0f, 0.0f)));
  EXPECT_EQ(1, std::count(hull.begin(), hull.end(), Point2f(0.0f, 1.0f)));
  EXPECT_EQ(1, std::count(hull.begin(), hull.end(), Point2f(1.0f, 1.0f)));
}

TEST(ConvexHull2dTest, TestEmptyHull) {
  // A set of colinear points should have an empty convex-hull.
  std::vector<Point2f> points;
  points.push_back(Point2f(0.0f, 0.0f));
  points.push_back(Point2f(2.0f, 2.0f));
  points.push_back(Point2f(4.0f, 4.0f));
  points.push_back(Point2f(6.0f, 6.0f));

  auto hull = ComputeConvexHull(points);
  EXPECT_EQ(0, hull.size());
}

}  // namespace geometry
}  // namespace seurat
