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

#include "seurat/tiler/subdivision.h"

#include <algorithm>
#include <array>
#include <random>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/util.h"
#include "seurat/geometry/fibonacci_sphere.h"
#include "seurat/tiler/point_set.h"
#include "seurat/tiler/tiler_test_util.h"

namespace seurat {
namespace tiler {
namespace {

using ion::math::Point3d;
using ion::math::Point3f;
using ion::math::Vector3f;

TEST(SubdivisionTest, TestPartitioning) {
  CubemapQuadtreeSubdivision pyramid(4);
  std::vector<Point3f> positions;

  // 50k points should be enough that at least one falls inside of each cell of
  // the quadtree.
  const int kPointCount = 50000;
  std::mt19937 random;
  std::uniform_real_distribution<float> radius(0.1f, 100.0f);
  for (int i = 0; i < kPointCount; ++i) {
    Point3d point = geometry::GenerateFibonacciSpherePoint(kPointCount, 0.0, i);
    positions.push_back(Point3f(point) * radius(random));
  }

  PointSet point_set;
  point_set.id = 0;
  point_set.positions = positions;

  pyramid.Init(point_set);

  ExpectWellFormedSubdivision(point_set, pyramid, 4);
}

TEST(SubdivisionTest, TestPartitioning_Hemisphere) {
  CubemapQuadtreeSubdivision pyramid(1);
  std::vector<Point3f> positions;

  const int kPointCount = 5000;
  std::mt19937 random;
  std::uniform_real_distribution<float> radius(0.1f, 100.0f);
  for (int i = 0; i < kPointCount; ++i) {
    // Generate points with positive-z only.
    Point3d point = geometry::GenerateFibonacciSpherePoint(kPointCount, 0.0, i);
    point[2] = std::fabs(point[2]);
    positions.push_back(Point3f(point) * radius(random));
  }

  PointSet point_set;
  point_set.id = 0;
  point_set.positions = positions;

  pyramid.Init(point_set);

  ExpectWellFormedSubdivision(point_set, pyramid, 4);
}

}  // namespace
}  // namespace tiler
}  // namespace seurat
