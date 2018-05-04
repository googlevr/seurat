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

#include "seurat/tiler/tiler_test_util.h"

#include <numeric>
#include <random>

#include "gtest/gtest.h"
#include "absl/types/span.h"
#include "seurat/base/util.h"
#include "seurat/geometry/mesh.h"
#include "seurat/geometry/mesh_util.h"
#include "seurat/geometry/plane.h"
#include "seurat/geometry/raytracer.h"

namespace seurat {
namespace tiler {

using geometry::Mesh;
using geometry::Plane3f;
using geometry::Raytracer;
using ion::math::Point3f;
using ion::math::Vector3f;

void FakeGeometrySolver::Init(const PointSet& point_set) {
  point_set_ = point_set;
}

void FakeGeometrySolver::InitializeModel(int point_index,
                                         GeometryModel* model) const {
  CHECK_GE(point_index, 0);
  CHECK_LT(point_index, point_set_.positions.size());
  Point3f center(point_index, 0.0f, 0.0f);
  Vector3f normal = Vector3f::AxisZ();
  model->center = center;
  model->normal = normal;
}

bool FakeGeometrySolver::FitModel(absl::Span<const int> point_indices,
                                  GeometryModel* model) const {
  std::vector<int> mutable_indices(point_indices.begin(), point_indices.end());
  std::vector<int>::iterator median =
      mutable_indices.begin() + (mutable_indices.size() / 2);
  std::nth_element(mutable_indices.begin(), median, mutable_indices.end());
  InitializeModel(*median, model);
  return true;
}

float FakeGeometrySolver::ComputeError(int point_index,
                                       const GeometryModel& model) const {
  CHECK_GE(point_index, 0);
  CHECK_LT(point_index, point_set_.positions.size());
  return std::fabs(point_index - model.center[0]);
}

float InvalidGeometrySolver::ComputeError(int point_index,
                                          const GeometryModel& model) const {
  return std::numeric_limits<float>::infinity();
}

void TestingPartitionerStage::Init(const PointSet& point_set) {
  point_set_ = point_set;
  init_count_++;
}

void TestingPartitionerStage::Run(const PointSet& point_set,
                                  absl::Span<BuildPartition> build_partitions) {
  CHECK_EQ(point_set_.id, point_set.id);
  CHECK_EQ(point_set_.positions.size(), point_set.positions.size());
  CHECK_EQ(point_set_.normals.size(), point_set.normals.size());
  CHECK_EQ(point_set_.colors.size(), point_set.colors.size());

  run_count_++;
}

void ExpectAllPointsPresent(const std::vector<BuildPartition>& partitioning,
                            int point_count) {
  std::vector<int> points_in_partitioning;
  for (const BuildPartition& bp : partitioning) {
    std::copy(bp.GetPointIndices().begin(), bp.GetPointIndices().end(),
              std::back_inserter(points_in_partitioning));
  }
  std::sort(points_in_partitioning.begin(), points_in_partitioning.end());

  std::vector<int> all_point_indices(point_count);
  std::iota(all_point_indices.begin(), all_point_indices.end(), 0);

  EXPECT_EQ(all_point_indices.size(), points_in_partitioning.size());
  EXPECT_EQ(all_point_indices, points_in_partitioning);
}

void ExpectNoDuplicatePoints(const std::vector<BuildPartition>& partitioning) {
  std::vector<int> point_indices;
  for (const auto& partition : partitioning) {
    std::copy(partition.GetPointIndices().begin(),
              partition.GetPointIndices().end(),
              std::back_inserter(point_indices));
  }
  std::sort(point_indices.begin(), point_indices.end());
  int num_unique = std::unique(point_indices.begin(), point_indices.end()) -
                   point_indices.begin();
  EXPECT_EQ(num_unique, point_indices.size());
}

float ComputeTotalError(const GeometrySolver& geometry_solver,
                        const std::vector<BuildPartition>& partitioning) {
  float total_error = 0.0f;
  for (const BuildPartition& bp : partitioning) {
    for (int point_index : bp.GetPointIndices()) {
      total_error += geometry_solver.ComputeError(point_index, bp.GetModel());
    }
  }
  return total_error;
}

std::vector<Point3f> GenerateUnitSpherePoints(int num_points) {
  std::mt19937 prng;
  std::normal_distribution<> normal_dist(0.0f, 1.0f);

  std::vector<Point3f> points;
  points.reserve(num_points);
  for (int i = 0; i < num_points; ++i) {
    // Generate a random point on the unit-sphere.
    Vector3f dir = Vector3f::Zero();
    while (!ion::math::Normalize(&dir)) {
      dir = Vector3f(normal_dist(prng), normal_dist(prng), normal_dist(prng));
    }
    points.push_back(Point3f::Zero() + dir);
  }
  return points;
}

void ExpectTilesCoverAllPoints(absl::Span<const Tile> tiles,
                               absl::Span<const ion::math::Point3f> points) {
  // Construct a raytracer for the tiles.
  Mesh mesh(0);
  for (const auto& tile : tiles) {
    geometry::AppendTriangleFan(tile.quad, {}, &mesh);
  }
  std::unique_ptr<Raytracer> raytracer = Raytracer::Build(mesh);

  // Verify that the ray from the origin to all points intersects a tile.
  for (const auto& point : points) {
    float t_hit;
    int triangle_hit_index;
    EXPECT_TRUE(raytracer->FindFirstHit(
        Point3f::Zero(), point - Point3f::Zero(), &t_hit, &triangle_hit_index));
  }
}

void ExpectTilesSphericalCap(Tiler* tiler, int point_set_id, int point_count,
                             const Vector3f& direction, float min_dot) {
  std::vector<Point3f> points = GenerateUnitSpherePoints(point_count);
  points.erase(std::remove_if(points.begin(), points.end(),
                              [=](const Point3f& point) {
                                float dot = ion::math::Dot(
                                    direction, point - Point3f::Zero());
                                return dot < min_dot;
                              }),
               points.end());

  const std::vector<Tile> tiles =
      tiler->Run(PointSet{point_set_id, points, {}, {}});

  ExpectTilesCoverAllPoints(tiles, points);
}

void ExpectTilesUnitSphere(Tiler* tiler, int point_set_id, int point_count) {
  return ExpectTilesSphericalCap(tiler, point_set_id, point_count,
                                 Vector3f::AxisZ(), -1.0f);
}

namespace {

void ExpectCounterClockwiseCellRailsContainAllPoints(
    const PointSet& all_points, const std::array<Vector3f, 4>& rails,
    absl::Span<const int> points_in_cell) {
  for (int i = 0; i < 4; ++i) {
    const Vector3f& v0 = rails[i];
    const Vector3f& v1 = rails[(i + 1) % 4];
    Vector3f inside = ion::math::Cross(v1, v0);
    Plane3f plane(Point3f::Zero(), inside);
    for (int point_index : points_in_cell) {
      const Point3f& point = all_points.positions[point_index];
      EXPECT_LE(0.0f, plane.SignedDistanceToPoint(point));
    }
  }
}

}  // namespace

void ExpectWellFormedSubdivision(const PointSet& points,
                                 const Subdivision& subdivision,
                                 int max_depth) {
  const int point_count = points.positions.size();

  std::vector<int> all_point_indices;
  std::vector<int> roots;
  subdivision.GetRoots(&roots);
  for (int cell : roots) {
    absl::Span<const int> point_indices = subdivision.GetPointsInCell(cell);
    all_point_indices.insert(all_point_indices.end(), point_indices.begin(),
                             point_indices.end());
  }
  // Verify that root cell contains all points.
  EXPECT_EQ(point_count, all_point_indices.size());

  // All point indices must exist exactly once.
  std::sort(all_point_indices.begin(), all_point_indices.end());
  for (int i = 0; i < point_count; ++i) {
    EXPECT_EQ(i, all_point_indices[i]);
  }

  // DFS over the subdivision, verifying that each child node contains a subset
  // of its parent's points.
  struct WorkItem {
    int cell;
    int depth;
  };
  std::vector<WorkItem> to_visit;
  to_visit.reserve(roots.size());
  for (int cell : roots) {
    to_visit.push_back({cell, 0});
  }
  std::vector<int> children;
  while (!to_visit.empty()) {
    WorkItem item = to_visit.back();
    int cur_node = item.cell;
    EXPECT_GE(max_depth, item.depth);
    to_visit.pop_back();

    children.clear();
    subdivision.GetChildren(cur_node, &children);

    ExpectCounterClockwiseCellRailsContainAllPoints(
        points, subdivision.GetRails(cur_node),
        subdivision.GetPointsInCell(cur_node));

    if (!children.empty()) {
      // The union over the points in each child node.
      std::vector<int> all_children_points;
      for (int child_node : children) {
        absl::Span<const int> child_node_points =
            subdivision.GetPointsInCell(child_node);
        all_children_points.insert(all_children_points.end(),
                                   child_node_points.begin(),
                                   child_node_points.end());
      }

      // Expect that the child nodes contain all the same points as the current
      // node.
      absl::Span<const int> cur_node_points =
          subdivision.GetPointsInCell(cur_node);
      std::vector<int> cur_node_points_sorted(cur_node_points.begin(),
                                              cur_node_points.end());
      std::sort(cur_node_points_sorted.begin(), cur_node_points_sorted.end());
      std::sort(all_children_points.begin(), all_children_points.end());
      EXPECT_EQ(cur_node_points_sorted, all_children_points);

      for (int child : children) {
        to_visit.push_back({child, item.depth + 1});
      }
    }
  }
}

}  // namespace tiler
}  // namespace seurat
