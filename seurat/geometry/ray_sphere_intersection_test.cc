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

#include "seurat/geometry/ray_sphere_intersection.h"

#include "ion/math/vector.h"
#include "gtest/gtest.h"

namespace seurat {
namespace geometry {
namespace {

using ion::math::Point3f;
using ion::math::Vector3f;

TEST(RaySphereIntersectionTest, TestIntersectRayStartingOutsideSphere) {
  const Point3f sphere_center = Point3f::Zero();
  float sphere_radius = 0.5f;
  const Point3f ray_start(0.0f, 0.0f, 1.0f);
  const Vector3f ray_direction(0.0f, 0.0f, -2.0f);

  float t_hit;
  EXPECT_TRUE(ComputeRaySphereIntersection(sphere_center, sphere_radius,
                                           ray_start, ray_direction, &t_hit));
  EXPECT_NEAR(0.25f, t_hit, 1e-5f);
}

TEST(RaySphereIntersectionTest, TestIntersectRayStartingInsideSphere) {
  const Point3f sphere_center = Point3f::Zero();
  float sphere_radius = 0.5f;
  const Point3f ray_start(0.0f, 0.0f, 0.0f);
  const Vector3f ray_direction(0.0f, 0.0f, 0.1f);

  float t_hit;
  EXPECT_TRUE(ComputeRaySphereIntersection(sphere_center, sphere_radius,
                                           ray_start, ray_direction, &t_hit));
  EXPECT_NEAR(5.0f, t_hit, 1e-5f);
}

TEST(RaySphereIntersectionTest, TestMissRayStartingOutsideSphere) {
  // Test with a ray which lies on a line intersecting the sphere.
  const Point3f sphere_center = Point3f::Zero();
  float sphere_radius = 0.5f;
  const Point3f ray_start(0.0f, 0.0f, 0.51);
  const Vector3f ray_direction(0.0f, 0.0f, 1.0f);

  float t_hit;
  EXPECT_FALSE(ComputeRaySphereIntersection(sphere_center, sphere_radius,
                                            ray_start, ray_direction, &t_hit));
}

TEST(RaySphereIntersectionTest, TestMissRay) {
  // Test with a ray which lies on a line which does not intersect the sphere.
  const Point3f sphere_center = Point3f::Zero();
  float sphere_radius = 0.5f;
  const Point3f ray_start(0.0f, 0.51f, 0.0f);
  const Vector3f ray_direction(0.0f, 0.0f, 1.0f);

  float t_hit;
  EXPECT_FALSE(ComputeRaySphereIntersection(sphere_center, sphere_radius,
                                            ray_start, ray_direction, &t_hit));
}

}  // namespace
}  // namespace geometry
}  // namespace seurat
