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

#include "seurat/tiler/rail_disk_solver.h"

#include <numeric>
#include <random>

#include "gtest/gtest.h"
#include "seurat/geometry/plane.h"
#include "seurat/tiler/point_set.h"

namespace seurat {
namespace tiler {
namespace {

using geometry::Plane3f;
using ion::math::Point3f;
using ion::math::Vector3f;
using seurat::base::Color3f;

double TotalError(const GeometrySolver& solver,
                  const std::vector<int>& point_indices,
                  const GeometryModel& model) {
  double total_err = 0.0;
  for (int point_index : point_indices) {
    total_err += solver.ComputeError(point_index, model);
  }
  return total_err;
}

int GetRootCellWithPoint(const Subdivision& subdivision, int point_index) {
  std::vector<int> roots;
  subdivision.GetRoots(&roots);
  int found_root = -1;
  for (int root : roots) {
    absl::Span<const int> points_in_root = subdivision.GetPointsInCell(root);
    if (std::find(points_in_root.begin(), points_in_root.end(), point_index) !=
        points_in_root.end()) {
      CHECK_EQ(found_root, -1)
          << "Malformed subdivision!  Point exists in multiple subtrees";
      found_root = root;
    }
  }

  CHECK_NE(found_root, -1) << "Point index not found in Subdivision";
  return found_root;
}

TEST(RailDiskSolverTest, InitializeSurfaceProxyAtEdge) {
  std::array<Point3f, 1> points = {{Point3f(10.0f, 0.0f, 0.0f)}};
  PointSet point_set{0, points, {}, {}};

  auto subdivision = std::make_shared<CubemapQuadtreeSubdivision>(3);
  subdivision->Init(point_set);
  int cell = GetRootCellWithPoint(*subdivision, 0);
  RailDiskSolver solver(0.0, subdivision, {1.0e-3f, 1.0e6f});
  solver.Init(point_set);

  GeometryModel model;
  model.cell = cell;
  solver.InitializeModel(0, &model);
  float error = solver.ComputeError(0, model);
  EXPECT_TRUE(std::isfinite(error));
  EXPECT_GT(std::numeric_limits<float>::max(), error);
  EXPECT_LE(0.0f, error);
}

TEST(RailDiskSolverTest, FitDisks_NotGrazingAngle) {
  const int kPointCount = 1000;

  Vector3f plane_vec1 = ion::math::Normalized(Vector3f(0.5f, 0.3f, 0.2f));
  Vector3f plane_vec2 = ion::math::Normalized(Vector3f(-4.0f, 4.9f, 0.3f));
  Vector3f plane_normal =
      ion::math::Normalized(ion::math::Cross(plane_vec1, plane_vec2));
  Point3f plane_pt = Point3f(0.0f, 3.0f, 4.0f) - plane_normal * 500.0f;

  std::vector<Point3f> points;
  std::vector<Vector3f> normals;
  std::vector<float> weights;
  // Generate points on the plane with random perturbations.
  std::mt19937 prng(0);
  // Standard normal distribution with mean = 0.0, sd = 1.0.
  std::normal_distribution<float> dist(0.0f, 1.0f);
  for (int i = 0; i < kPointCount; ++i) {
    const Point3f position = plane_pt + plane_vec1 * dist(prng) * 50.0f +
                             plane_vec2 * dist(prng) * 50.0f +
                             plane_normal * dist(prng);
    const Vector3f normal = ion::math::Normalized(
        plane_normal + Vector3f(dist(prng), dist(prng), dist(prng)) * 0.25f);
    points.push_back(position);
    normals.push_back(normal);
    weights.push_back(std::fabs(dist(prng) + 1.0f));
  }

  PointSet point_set;
  point_set.id = 0;
  point_set.positions = points;
  point_set.normals = normals;
  point_set.weights = weights;

  auto subdivision = std::make_shared<CubemapQuadtreeSubdivision>(3);
  subdivision->Init(point_set);
  // Use the highest-level (largest) cell which is not the root node.
  int cell = GetRootCellWithPoint(*subdivision, 0);

  RailDiskSolver solver(0.01f, subdivision, {1.0e-3f, 1.0e8f});
  solver.Init(point_set);

  std::vector<int> all_point_indices(points.size());
  std::iota(all_point_indices.begin(), all_point_indices.end(), 0);

  // Note that the disk-fitting solver may have multiple extrema, so choice of
  // the initial model estimate does matter.
  //
  // Starting with a perturbed version of the underlying plane's normal is good
  // enough to find a local minimum which is adequate for testing purposes.
  GeometryModel initial_model_estimate;
  initial_model_estimate.cell = cell;
  initial_model_estimate.center = plane_pt;
  initial_model_estimate.normal =
      ion::math::Normalized(plane_normal + Vector3f(0.1f, -0.05f, -0.1f));
  GeometryModel fit_model = initial_model_estimate;
  bool success = solver.FitModel(all_point_indices, &fit_model);
  EXPECT_TRUE(success);
  Plane3f fit_plane = fit_model.GetPlane();

  if (ion::math::Dot(plane_normal, fit_plane.GetNormal()) < 0.0f) {
    plane_normal *= -1.0f;
  }

  EXPECT_LT(0.9f, ion::math::Dot(plane_normal, fit_plane.GetNormal()));

  float err_at_fit_model = TotalError(solver, all_point_indices, fit_model);

  // At the very least, the resulting GeometryModel should be near locally
  // optimal (guaranteeing a global optimum is intractable).
  //
  // So, we can expect that the total error for perturbed versions of the fit
  // GeometryModel should be greater.

  for (int i = 0; i < 100; ++i) {
    Vector3f random_direction =
        ion::math::Normalized(Vector3f{dist(prng), dist(prng), dist(prng)});
    GeometryModel translated_model = fit_model;
    translated_model.center += random_direction * 1.0f;
    float err_at_perturbed =
        TotalError(solver, all_point_indices, translated_model);
    EXPECT_GE(err_at_perturbed, err_at_fit_model);
  }

  for (int i = 0; i < 100; ++i) {
    Vector3f random_direction =
        ion::math::Normalized(Vector3f{dist(prng), dist(prng), dist(prng)});
    GeometryModel perturbed_model = fit_model;
    perturbed_model.normal =
        ion::math::Normalized(fit_model.normal + random_direction * 0.1f);
    float err_at_perturbed =
        TotalError(solver, all_point_indices, perturbed_model);
    EXPECT_GE(err_at_perturbed, err_at_fit_model);
  }
}

void TestGrazingAngleQuad(const std::array<Point3f, 4>& quad) {
  // Generate an 11x11 grid of points on the quad.
  std::vector<Point3f> points;
  for (int u = 0; u <= 10; ++u) {
    const Point3f interp0 = (quad[1] - quad[0]) * (u / 10.0f) + quad[0];
    const Point3f interp1 = (quad[2] - quad[3]) * (u / 10.0f) + quad[2];
    for (int v = 0; v <= 10; ++v) {
      Point3f point = (interp1 - interp0) * (v / 10.0f) + interp0;
      points.push_back(point);
    }
  }
  // Index into 'points' of a point in the middle of the quad.
  const int middle_point_index = 11 * 6 + 6;

  PointSet point_set;
  point_set.id = 0;
  point_set.positions = points;

  auto subdivision = std::make_shared<CubemapQuadtreeSubdivision>(3);
  subdivision->Init(point_set);
  // Use the highest-level (largest) cell which is not the root node.
  //
  // This should correspond to a cube face.
  int cell = GetRootCellWithPoint(*subdivision, middle_point_index);

  RailDiskSolver solver(0.01f, subdivision, {1.0e-3f, 1.0e8f});
  solver.Init(point_set);

  std::vector<int> all_point_indices(points.size());
  std::iota(all_point_indices.begin(), all_point_indices.end(), 0);

  GeometryModel model;
  model.cell = cell;
  // Initialize the model with a point in the middle.
  solver.InitializeModel(11 * 6 + 6, &model);
  bool success = solver.FitModel(all_point_indices, &model);
  EXPECT_TRUE(success);

  // The resulting model's plane should intersect the cell's rails.
  Plane3f plane = model.GetPlane();
  std::array<Vector3f, 4> rails = subdivision->GetRails(cell);
  for (const auto& rail : rails) {
    float t_hit;
    EXPECT_TRUE(plane.IntersectRay(Point3f::Zero(), rail, &t_hit));
  }
}

TEST(RailDiskSolverTest, FitDisks_GrazingAngle) {
  // A quad at a grazing angle within the -z cube face.
  std::array<Point3f, 4> quad = {
      {Point3f(-5.0f, -0.1f, -10.0f), Point3f(5.0f, -0.1f, -10.0f),
       Point3f(5.0f, -0.2f, -20.0f), Point3f(-5.0f, -0.2f, -20.0f)}};
  TestGrazingAngleQuad(quad);
}

TEST(RailDiskSolverTest, FitDisks_Perpendicular) {
  // A quad parallel to the xz plane.
  std::array<Point3f, 4> quad = {
      {Point3f(-5.0f, 0.0f, -10.0f), Point3f(5.0f, 0.0f, -10.0f),
       Point3f(5.0f, 0.0f, -20.0f), Point3f(-5.0f, 0.0f, -20.0f)}};
  TestGrazingAngleQuad(quad);
}

}  // namespace
}  // namespace tiler
}  // namespace seurat
