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

#include "seurat/geometry/kdtree.h"

#include <vector>

#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace seurat {
namespace geometry {

using ion::math::Point2f;
using ion::math::Point3f;

namespace {

std::vector<Point3f> BuildUnitSquarePointCloud() {
  std::vector<Point3f> points;

  // Create a point cloud consisting of vertices of a 2x2 square in the
  // xy-plane.
  points.push_back({0.0f, 0.0f, 0.0f});
  points.push_back({0.0f, 2.0f, 0.0f});
  points.push_back({2.0f, 0.0f, 0.0f});
  points.push_back({2.0f, 2.0f, 0.0f});
  return points;
}

}  // namespace

TEST(KdTreeTest, KnnSearch) {
  auto points = BuildUnitSquarePointCloud();
  const KdTree<3> kdtree(points);
  std::vector<int> result;

  result.clear();
  kdtree.KnnSearch(points[0], 1, &result);
  EXPECT_EQ((std::vector<int>{0}), result);

  result.clear();
  kdtree.KnnSearch(points[1], 1, &result);
  EXPECT_EQ((std::vector<int>{1}), result);

  result.clear();
  kdtree.KnnSearch(points[2], 1, &result);
  EXPECT_EQ((std::vector<int>{2}), result);

  result.clear();
  kdtree.KnnSearch(points[3], 1, &result);
  EXPECT_EQ((std::vector<int>{3}), result);

  result.clear();
  kdtree.KnnSearch(points[0], 2, &result);
  EXPECT_THAT(result,
              ::testing::UnorderedElementsAreArray(std::vector<int>{0, 1}));

  result.clear();
  kdtree.KnnSearch(points[1], 2, &result);
  EXPECT_THAT(result,
              ::testing::UnorderedElementsAreArray(std::vector<int>{1, 0}));

  result.clear();
  kdtree.KnnSearch(points[2], 2, &result);
  EXPECT_THAT(result,
              ::testing::UnorderedElementsAreArray(std::vector<int>{2, 0}));

  result.clear();
  kdtree.KnnSearch(points[3], 2, &result);
  EXPECT_THAT(result,
              ::testing::UnorderedElementsAreArray(std::vector<int>{3, 1}));
}

TEST(KdTreeTest, NNSearch) {
  auto points = BuildUnitSquarePointCloud();
  const KdTree<3> kdtree(points);

  int result;
  kdtree.NNSearch(points[0], &result);
  EXPECT_EQ(0, result);

  kdtree.NNSearch(points[1], &result);
  EXPECT_EQ(1, result);

  kdtree.NNSearch(points[2], &result);
  EXPECT_EQ(2, result);

  kdtree.NNSearch(points[3], &result);
  EXPECT_EQ(3, result);
}

TEST(KdTreeTest, KnnSearch_Empty) {
  std::vector<Point3f> no_points;
  const KdTree<3> kdtree(no_points);
  std::vector<int> result;

  result.clear();
  kdtree.KnnSearch(Point3f::Zero(), 1, &result);
  EXPECT_TRUE(result.empty());
}

TEST(KdTreeTest, KnnSearch_TooFewPoints) {
  auto points = BuildUnitSquarePointCloud();
  const KdTree<3> kdtree(points);
  std::vector<int> result;

  result.clear();
  kdtree.KnnSearch(Point3f::Zero(), points.size() * 3, &result);
  EXPECT_EQ(points.size(), result.size());
}

TEST(KdTreeTest, RadiusSearch) {
  auto points = BuildUnitSquarePointCloud();
  const KdTree<3> kdtree(points);
  std::vector<int> result;

  result.clear();
  kdtree.RadiusSearch(points[0], 3.9f, false, &result);
  EXPECT_EQ((std::vector<int>{0}), result);

  result.clear();
  kdtree.RadiusSearch(points[1], 3.9f, false, &result);
  EXPECT_EQ((std::vector<int>{1}), result);

  result.clear();
  kdtree.RadiusSearch(points[2], 3.9f, false, &result);
  EXPECT_EQ((std::vector<int>{2}), result);

  result.clear();
  kdtree.RadiusSearch(points[3], 3.9f, false, &result);
  EXPECT_EQ((std::vector<int>{3}), result);

  result.clear();
  kdtree.RadiusSearch(points[0], 4.1f, false, &result);
  EXPECT_THAT(result,
              ::testing::UnorderedElementsAreArray(std::vector<int>{0, 1, 2}));

  result.clear();
  kdtree.RadiusSearch(points[1], 4.1f, false, &result);
  EXPECT_THAT(result,
              ::testing::UnorderedElementsAreArray(std::vector<int>{0, 1, 3}));

  result.clear();
  kdtree.RadiusSearch(points[2], 4.1f, false, &result);
  EXPECT_THAT(result,
              ::testing::UnorderedElementsAreArray(std::vector<int>{0, 2, 3}));

  result.clear();
  kdtree.RadiusSearch(points[3], 4.1f, false, &result);
  EXPECT_THAT(result,
              ::testing::UnorderedElementsAreArray(std::vector<int>{1, 2, 3}));
}

TEST(KdTreeTest, Test_2D_KdTree) {
  std::vector<Point2f> points;
  // Create a point cloud consisting of vertices of a 2x2 square.
  points.push_back({0.0f, 0.0f});
  points.push_back({0.0f, 2.0f});
  points.push_back({2.0f, 0.0f});
  points.push_back({2.0f, 2.0f});

  const KdTree<2> kdtree(points);

  std::vector<int> result;
  kdtree.KnnSearch(points[1], 2, &result);
  EXPECT_THAT(result,
              ::testing::UnorderedElementsAreArray(std::vector<int>{1, 0}));

  result.clear();

  kdtree.RadiusSearch(points[2], 4.1f, false, &result);
  EXPECT_THAT(result,
              ::testing::UnorderedElementsAreArray(std::vector<int>{0, 2, 3}));
}

}  // namespace geometry
}  // namespace seurat
