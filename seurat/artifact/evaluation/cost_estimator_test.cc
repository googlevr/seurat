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

#include "seurat/artifact/evaluation/cost_estimator.h"

#include "ion/math/matrix.h"
#include "ion/math/matrixutils.h"
#include "ion/math/range.h"
#include "ion/math/rangeutils.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/geometry/mesh.h"
#include "seurat/geometry/mesh_util.h"
#include "seurat/mesh/mesh_component_util.h"

namespace seurat {
namespace artifact {

using geometry::Mesh;
using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Range3f;
using ion::math::Vector3f;

namespace {

// Returns a mesh consisting of an origin-centered axis-aligned box with the
// side-lengths specified by |diameter|.
Mesh GetCubeMesh(const Vector3f& diameter) {
  Mesh cube =
      mesh::MeshComponentUtil::ToMesh(*mesh::MeshComponentUtil::CreateCube({}));
  cube.TransformPositions(ion::math::ScaleMatrixH(diameter));
  return cube;
}

}  // namespace

TEST(CostEstimatorTest, TestEmpty) {
  CostEstimator cost_estimator = CostEstimator::Build({});
  EXPECT_EQ(0.0f, cost_estimator.EstimateCost(
                      Range3f::BuildWithSize({-1.0f, -1.0f, -1.0f},
                                             {2.0f, 2.0f, 2.0f}),
                      3, 3));
}

TEST(CostEstimatorTest, TestSingleLayerOpaqueGeometry) {
  Mesh mesh = GetCubeMesh({10.0f, 10.0f, 10.0f});
  EXPECT_EQ(1, mesh.GetTextureCount());

  CostEstimator cost_estimator = CostEstimator::Build({{mesh, 2.4f, true}});

  EXPECT_NEAR(
      2.4f,
      cost_estimator.EstimateCost(
          Range3f::BuildWithSize({-1.0f, -1.0f, -1.0f}, {2.0f, 2.0f, 2.0f}), 3,
          3),
      1e-5f);
}

TEST(CostEstimatorTest, TestMultipleLayerOpaqueGeometry) {
  Mesh mesh = GetCubeMesh({10.0f, 10.0f, 10.0f});
  mesh.AppendMesh(GetCubeMesh({11.0f, 9.0f, 10.0f}));

  CostEstimator cost_estimator = CostEstimator::Build({{mesh, 2.4f, true}});

  EXPECT_NEAR(
      2.4f,
      cost_estimator.EstimateCost(
          Range3f::BuildWithSize({-1.0f, -1.0f, -1.0f}, {2.0f, 2.0f, 2.0f}), 3,
          3),
      1e-5f);
}

TEST(CostEstimatorTest, TestMultipleMeshOpaqueGeometry) {
  Mesh mesh1 = GetCubeMesh({10.0f, 10.0f, 10.0f});
  Mesh mesh2 = GetCubeMesh({15.0f, 15.0f, 15.0f});

  CostEstimator cost_estimator = CostEstimator::Build({
      {mesh1, 2.4f, true},
      {mesh2, 3.0f, true},
  });

  // The outer mesh, mesh2, should have no impact on total-cost, since all rays
  // must intersect mesh1 first.
  EXPECT_NEAR(
      2.4f,
      cost_estimator.EstimateCost(
          Range3f::BuildWithSize({-1.0f, -1.0f, -1.0f}, {2.0f, 2.0f, 2.0f}), 3,
          3),
      1e-5f);
}

TEST(CostEstimatorTest, TestSingleLayerTransparentGeometry) {
  Mesh mesh = GetCubeMesh({10.0f, 10.0f, 10.0f});

  CostEstimator cost_estimator = CostEstimator::Build({
      {mesh, 2.4f, false},
  });

  EXPECT_NEAR(
      2.4f,
      cost_estimator.EstimateCost(
          Range3f::BuildWithSize({-1.0f, -1.0f, -1.0f}, {2.0f, 2.0f, 2.0f}), 3,
          3),
      1e-5f);
}

TEST(CostEstimatorTest, TestMultipleLayerTransparentGeometry) {
  Mesh mesh = GetCubeMesh({10.0f, 10.0f, 10.0f});
  mesh.AppendMesh(GetCubeMesh({11.0f, 11.0f, 11.0f}));

  CostEstimator cost_estimator = CostEstimator::Build({
      {mesh, 2.4f, false},
  });

  EXPECT_NEAR(
      2.4f + 2.4f,
      cost_estimator.EstimateCost(
          Range3f::BuildWithSize({-1.0f, -1.0f, -1.0f}, {2.0f, 2.0f, 2.0f}), 3,
          3),
      1e-5f);
}

TEST(CostEstimatorTest, TestMultipleMeshTransparentGeometry) {
  // Two cuboids with z-fighting along the shared xy-parallel sides.
  Mesh mesh1 = GetCubeMesh({10.0f, 10.0f, 10.0f});
  Mesh mesh2 = GetCubeMesh({15.0f, 5.0f, 10.0f});

  CostEstimator cost_estimator = CostEstimator::Build({
      {mesh1, 2.4f, false},
      {mesh2, 1.3f, false},
  });

  float expected_average_cost = 2.4f + 1.3f;

  EXPECT_NEAR(
      expected_average_cost,
      cost_estimator.EstimateCost(
          Range3f::BuildWithSize({-1.0f, -1.0f, -1.0f}, {2.0f, 2.0f, 2.0f}), 3,
          3),
      1e-5f);
}

TEST(CostEstimatorTest, TestMultipleMeshMixedGeometry) {
  const float kBaselineSize = 10000.0f;
  Mesh opaque_mesh = GetCubeMesh({kBaselineSize, kBaselineSize, kBaselineSize});

  for (int axis_to_test = 0; axis_to_test < 3; ++axis_to_test) {
    // The transparent mesh should be visible only on 4/6 sides.
    //
    // The axis_to_test determines which sides of this mesh will be visible.
    Vector3f transparent_mesh_diameter(
        kBaselineSize - 1.0f, kBaselineSize - 1.0f, kBaselineSize - 1.0f);
    transparent_mesh_diameter[axis_to_test] = kBaselineSize + 1.0f;

    Mesh transparent_mesh = GetCubeMesh(transparent_mesh_diameter);

    CostEstimator cost_estimator = CostEstimator::Build({
        {opaque_mesh, 2.4f, true},
        {transparent_mesh, 1.3f, false},
    });

    float expected_average_cost = 2.4f + 1.3f * 4.0f / 6.0f;

    EXPECT_NEAR(
        expected_average_cost,
        cost_estimator.EstimateCost(
            Range3f::BuildWithSize({-1.0f, -1.0f, -1.0f}, {2.0f, 2.0f, 2.0f}),
            1, 6000),
        1e-3f);
  }
}

TEST(CostEstimatorTest, TestMultipleMeshMixedGeometry_EstimateRayCost) {
  Mesh opaque_mesh = GetCubeMesh({1000.0f, 1000.0f, 1000.0f});

  // The transparent mesh should be visible only on 4/6 sides.
  Mesh transparent_mesh = GetCubeMesh({1001.0f, 999.0f, 999.0f});

  CostEstimator cost_estimator = CostEstimator::Build({
      {opaque_mesh, 2.4f, true},
      {transparent_mesh, 1.3f, false},
  });

  // The underlying implementation uses Embree to trace rays, which may report
  // duplicate intersections for rays which hit the center of a quad where its
  // two triangles meet.  So, use a slightly-perturbed ray origin to avoid that
  // edge-case.
  Point3f ray_origin(0.1f, 0.2f, 0.3f);

  EXPECT_EQ(2.4f,
            cost_estimator.EstimateRayCost(ray_origin, Vector3f::AxisX()));
  EXPECT_EQ(2.4f,
            cost_estimator.EstimateRayCost(ray_origin, -Vector3f::AxisX()));
  EXPECT_EQ(2.4f + 1.3f,
            cost_estimator.EstimateRayCost(ray_origin, Vector3f::AxisY()));
  EXPECT_EQ(2.4f + 1.3f,
            cost_estimator.EstimateRayCost(ray_origin, -Vector3f::AxisY()));
  EXPECT_EQ(2.4f + 1.3f,
            cost_estimator.EstimateRayCost(ray_origin, Vector3f::AxisZ()));
  EXPECT_EQ(2.4f + 1.3f,
            cost_estimator.EstimateRayCost(ray_origin, -Vector3f::AxisZ()));
}

}  // namespace artifact
}  // namespace seurat
