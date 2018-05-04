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

#include "seurat/geometry/fibonacci_sphere.h"

#include <cmath>
#include <random>

#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "gtest/gtest.h"

namespace seurat {
namespace geometry {
namespace {

using ion::math::Point3d;
using ion::math::Vector3d;

int CountPointsInHalfspace(const std::vector<Point3d>& points,
                           const Vector3d& halfspace) {
  int num_points = 0;
  for (const auto& pt : points) {
    if (ion::math::Dot(pt - Point3d::Zero(), halfspace) >= 0.0f) {
      num_points++;
    }
  }
  return num_points;
}

TEST(FibonacciSphereTest, GenerateFibonacciSpherePoint) {
  const int kNumPoints = 500;
  const double kScrambler = 1.0 * M_PI;
  std::vector<Point3d> points(kNumPoints);

  for (int i = 0; i < kNumPoints; ++i) {
    points[i] = GenerateFibonacciSpherePoint(kNumPoints, kScrambler, i);

    EXPECT_NEAR(1.0f, ion::math::Length(points[i] - Point3d::Zero()), 1.0e-4f);
  }

  // Verify that the number of points in each hemisphere is near kNumPoints / 2.
  //
  // Note that there is a lot of slack in these tests, since the resulting
  // point set is not expected to be perfectly uniform.
  const float kHalf = kNumPoints / 2.0f;
  const float kSlack = kNumPoints * 1.0e-2f;

  EXPECT_NEAR(kHalf, CountPointsInHalfspace(points, Vector3d::AxisX()), kSlack);
  EXPECT_NEAR(kHalf, CountPointsInHalfspace(points, Vector3d::AxisY()), kSlack);
  EXPECT_NEAR(kHalf, CountPointsInHalfspace(points, Vector3d::AxisZ()), kSlack);
  EXPECT_NEAR(kHalf, CountPointsInHalfspace(points, -Vector3d::AxisX()),
              kSlack);
  EXPECT_NEAR(kHalf, CountPointsInHalfspace(points, -Vector3d::AxisY()),
              kSlack);
  EXPECT_NEAR(kHalf, CountPointsInHalfspace(points, -Vector3d::AxisZ()),
              kSlack);
}

TEST(FibonacciSphereTest, MillionsOfPoints) {
  // 10 million points.
  const int kNumPoints = 10000000;
  const double kScrambler = 1.0 * M_PI;
  std::vector<Point3d> points(kNumPoints);

  for (int i = 0; i < kNumPoints; ++i) {
    points[i] = GenerateFibonacciSpherePoint(kNumPoints, kScrambler, i);

    EXPECT_NEAR(1.0f, ion::math::Length(points[i] - Point3d::Zero()), 1.0e-4f);
  }
}

TEST(FibonacciSphereTest, InverseMapping) {
  const int kNumPoints = 1000;
  const float kScrambler = 0.0f;

  for (int i = 0; i < kNumPoints; ++i) {
    Point3d point = GenerateFibonacciSpherePoint(kNumPoints, kScrambler, i);
    int inverse = InverseFibonacciSphereMapping(
        kNumPoints, kScrambler, Point3d(point) - Point3d::Zero());
    EXPECT_EQ(i, inverse);
  }
}

TEST(FibonacciSphereTest, InverseMappingWithScrambler) {
  const int kNumPoints = 1000;

  std::mt19937 random;
  std::uniform_real_distribution<double> angle(0, 2 * M_PI);
  for (int i = 0; i < kNumPoints; ++i) {
    double kScrambler = angle(random);

    Point3d point = GenerateFibonacciSpherePoint(kNumPoints, kScrambler, i);
    int inverse = InverseFibonacciSphereMapping(
        kNumPoints, kScrambler, Point3d(point) - Point3d::Zero());
    EXPECT_EQ(i, inverse);
  }
}

}  // namespace
}  // namespace geometry
}  // namespace seurat
