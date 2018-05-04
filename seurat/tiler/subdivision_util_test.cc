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

#include "seurat/tiler/subdivision_util.h"

#include <algorithm>
#include <array>
#include <random>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "absl/types/span.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/base/util.h"
#include "seurat/geometry/fibonacci_sphere.h"
#include "seurat/geometry/plane.h"
#include "seurat/tiler/point_set.h"
#include "seurat/tiler/subdivision.h"

namespace seurat {
namespace tiler {
namespace {

using ion::math::Point3f;

TEST(SubdivisionUtilTest, GetCellsInDepthRange) {
  CubemapQuadtreeSubdivision pyramid(5);
  // Initialize with unit sphere points.
  std::vector<Point3f> positions;
  const int kPointCount = 5000;
  for (int i = 0; i < kPointCount; ++i) {
    Point3f point(geometry::GenerateFibonacciSpherePoint(kPointCount, 0.0, i));
    positions.push_back(point);
  }

  PointSet point_set;
  point_set.id = 0;
  point_set.positions = positions;
  pyramid.Init(point_set);

  std::vector<int> cells;
  GetCellsInDepthRange(pyramid, 0, 0, &cells);
  EXPECT_FALSE(cells.empty());

  std::vector<int> roots;
  pyramid.GetRoots(&roots);
  std::sort(roots.begin(), roots.end());
  std::sort(cells.begin(), cells.end());
  EXPECT_EQ(cells, roots);

  std::vector<int> level_1;
  for (int root : roots) {
    pyramid.GetChildren(root, &level_1);
  }
  GetCellsInDepthRange(pyramid, 1, 1, &cells);
  std::sort(level_1.begin(), level_1.end());
  std::sort(cells.begin(), cells.end());
  EXPECT_EQ(level_1, cells);

  std::vector<int> roots_and_level_1;
  roots_and_level_1.insert(roots_and_level_1.end(), roots.begin(), roots.end());
  roots_and_level_1.insert(roots_and_level_1.end(), level_1.begin(),
                           level_1.end());
  GetCellsInDepthRange(pyramid, 0, 1, &cells);
  std::sort(roots_and_level_1.begin(), roots_and_level_1.end());
  std::sort(cells.begin(), cells.end());
  EXPECT_EQ(roots_and_level_1, cells);
}

}  // namespace
}  // namespace tiler
}  // namespace seurat
