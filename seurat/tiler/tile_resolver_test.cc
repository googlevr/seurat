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

#include "seurat/tiler/tile_resolver.h"

#include <algorithm>
#include <random>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/geometry/fibonacci_sphere.h"
#include "seurat/geometry/triangle.h"
#include "seurat/tiler/build_partition.h"
#include "seurat/tiler/point_set.h"
#include "seurat/tiler/tile.h"

namespace seurat {
namespace tiler {
namespace {

using geometry::Plane3f;
using geometry::Triangle3f;
using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Vector3f;

void ValidateTile(const Tile& tile, const PointSet& points,
                  const BuildPartition& partition) {
  const Plane3f& plane = partition.GetModel().GetPlane();

  // Verify that all points are within the cone from the origin to the tile's
  // quad.
  //
  // This guarantees that the projection of these points onto the tile's plane
  // will be within the bounds of the tile.
  //
  // Note that this is deliberately agnostic to the tile's winding order.
  std::vector<Plane3f> halfspaces;
  for (int i = 0; i < 4; ++i) {
    const Vector3f& v0 = tile.quad[i] - Point3f::Zero();
    const Vector3f& v1 = tile.quad[(i + 1) % 4] - Point3f::Zero();
    const Point3f& p2 = tile.quad[(i + 2) % 4];
    Vector3f inside = ion::math::Cross(v1, v0);
    Plane3f plane(Point3f::Zero(), inside);
    if (plane.SignedDistanceToPoint(p2) < 0.0f) {
      plane = plane.GetReversePlane();
    }
    halfspaces.push_back(plane);
  }
  for (int point_index : partition.GetPointIndices()) {
    const Point3f& point = points.positions[point_index];
    for (const auto& halfspace : halfspaces) {
      // Add some slack for points on the edge of the quad.
      EXPECT_LE(-1e-3f, halfspace.SignedDistanceToPoint(point));
    }
  }
  // Verify the quad is on the plane.
  EXPECT_NEAR(0.0f, plane.SignedDistanceToPoint(tile.quad[0]), 1e-3f);
  EXPECT_NEAR(0.0f, plane.SignedDistanceToPoint(tile.quad[1]), 1e-3f);
  EXPECT_NEAR(0.0f, plane.SignedDistanceToPoint(tile.quad[2]), 1e-3f);
  EXPECT_NEAR(0.0f, plane.SignedDistanceToPoint(tile.quad[3]), 1e-3f);
}

void TestTileResolverWithEasyTestCase(TileResolver* resolver) {
  std::vector<Point3f> points;
  std::mt19937 random;
  std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
  Point3f offset(1.0f, 10.0f, -2.0f);

  for (int i = 0; i < 20; ++i) {
    Point3f random_point(dist(random), dist(random), dist(random));
    random_point += offset;
    points.push_back(random_point);
  }

  BuildPartition partition;
  for (int i = 0; i < points.size(); ++i) {
    partition.AddPoint(i, 0.0f /* unused */);
  }
  partition.GetMutableModel()->center = offset;
  partition.GetMutableModel()->normal = Vector3f::AxisY();

  PointSet point_set;
  point_set.positions = points;
  resolver->Init(point_set);

  Tile tile;
  EXPECT_TRUE(resolver->Resolve(partition, &tile));
  ValidateTile(tile, point_set, partition);
}

void TestTileResolverWithDegeneratePlane(TileResolver* resolver) {
  std::vector<Point3f> points;
  std::mt19937 random;
  std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
  Point3f offset(1.0f, 10.0f, -2.0f);

  for (int i = 0; i < 20; ++i) {
    Point3f random_point(dist(random), dist(random), dist(random));
    random_point += offset;
    points.push_back(random_point);
  }

  // Verify failure if the partition plane intersects the origin.
  BuildPartition bad_partition;
  for (int i = 0; i < points.size(); ++i) {
    bad_partition.AddPoint(i, 0.0f /* unused */);
  }
  // This GeometryModel intersects the origin.
  bad_partition.GetMutableModel()->center = Point3f::Zero();
  bad_partition.GetMutableModel()->normal = Vector3f::AxisY();

  PointSet point_set;
  point_set.positions = points;
  resolver->Init(point_set);

  Tile tile;
  EXPECT_FALSE(resolver->Resolve(bad_partition, &tile));
}

TEST(TileResolverTest, RailTileResolver) {
  auto subdivision = std::make_shared<CubemapQuadtreeSubdivision>(3);
  RailTileResolver resolver(subdivision);
  const int kCell = 9;

  std::vector<Point3f> points;

  const float kRadius = 10.0f;

  const int kNumPoints = 1000;
  for (int i = 0; i < kNumPoints; ++i) {
    Point3f point(geometry::GenerateFibonacciSpherePoint(kNumPoints, 0.0, i));
    points.push_back(point * kRadius);
  }

  PointSet point_set;
  point_set.id = 0;
  point_set.positions = points;
  subdivision->Init(point_set);

  BuildPartition partition;
  for (int point_index : subdivision->GetPointsInCell(kCell)) {
    partition.AddPoint(point_index, 0.0f /* unused */);
  }
  partition.GetMutableModel()->cell = kCell;
  partition.GetMutableModel()->center =
      points[partition.GetPointIndices().front()];
  partition.GetMutableModel()->normal = ion::math::Normalized(
      partition.GetMutableModel()->center - Point3f::Zero());

  resolver.Init(point_set);

  Tile tile;
  EXPECT_TRUE(resolver.Resolve(partition, &tile));
  ValidateTile(tile, point_set, partition);
}

}  // namespace
}  // namespace tiler
}  // namespace seurat
