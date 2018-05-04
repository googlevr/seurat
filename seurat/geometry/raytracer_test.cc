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

#include "seurat/geometry/raytracer.h"

#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "gtest/gtest.h"
#include "seurat/geometry/mesh.h"
#include "seurat/geometry/mesh_util.h"

namespace seurat {
namespace geometry {
namespace {

using ion::math::Point3f;
using ion::math::Vector3f;

TEST(RaytracerTest, TestIntersectRayQuad) {
  // A unit-sized square in the z=-1 plane.
  const std::vector<Point3f> vertex_array = {
      Point3f(0.0f, 0.0f, -1.0f),
      Point3f(1.0f, 0.0f, -1.0f),
      Point3f(1.0f, 1.0f, -1.0f),
      Point3f(0.0f, 1.0f, -1.0f),
  };
  const std::vector<int> index_buffer = {0, 1, 2, 0, 2, 3};
  std::unique_ptr<Raytracer> raytracer =
      Raytracer::Build(vertex_array, index_buffer);

  float t_hit;
  int triangle_index;
  bool hit = raytracer->FindFirstHit(Point3f(0.5f, 0.5f, 0.5f),
                                     Vector3f(0.0f, 0.0f, -1.0f), &t_hit,
                                     &triangle_index);
  EXPECT_TRUE(hit);
  EXPECT_GT(1.0e-5f, std::fabs(t_hit - 1.5f));
}

TEST(RaytracerTest, TestNoIntersection) {
  // A unit-sized square in the z=-1 plane.
  const std::vector<Point3f> vertex_array = {
      Point3f(0.0f, 0.0f, -1.0f),
      Point3f(1.0f, 0.0f, -1.0f),
      Point3f(1.0f, 1.0f, -1.0f),
      Point3f(0.0f, 1.0f, -1.0f),
  };
  const std::vector<int> index_buffer = {0, 1, 2, 0, 2, 3};
  std::unique_ptr<Raytracer> raytracer =
      Raytracer::Build(vertex_array, index_buffer);

  float t_hit;
  int triangle_index;
  bool hit = raytracer->FindFirstHit(Point3f(0.5f, 0.5f, -3.0f),
                                     Vector3f(0.0f, 0.0f, -1.0f), &t_hit,
                                     &triangle_index);
  EXPECT_FALSE(hit);
}

// Builds a scene consisting of 3 overlapping quads:
//  1. A square at z=0.
//  2. A rectangle at z=0.1 covering the right half (x > 0.5) of the square.
//  3. A rectangle at z=0.2 on the top half (y > 0.5).
//
// The following diagram indicates the number of layers of geometry at different
// locations.
//
// ^ +y
// |
//
// +----+----+
// |  2 |  3 |
// |    |    |
// +----+----+
// |  1 |  2 |
// |    |    |
// +---------+ ->
//             +x
Mesh MakeTestScene() {
  Mesh mesh(0);
  // The unit-square in the XY plane.
  AppendTriangleFan(
      {
          Point3f(0.0f, 0.0f, 0.0f),  //
          Point3f(0.0f, 1.0f, 0.0f),  //
          Point3f(1.0f, 1.0f, 0.0f),  //
          Point3f(1.0f, 0.0f, 0.0f)   //
      },
      {}, &mesh);

  // The unit-square, clipped to x > 0.5 and translated to z = 0.1.
  AppendTriangleFan(
      {
          Point3f(0.5f, 0.0f, 0.1f),  //
          Point3f(0.5f, 1.0f, 0.1f),  //
          Point3f(1.0f, 1.0f, 0.1f),  //
          Point3f(1.0f, 0.0f, 0.1f)   //
      },
      {}, &mesh);

  // The unit-square, clipped to y > 0.5 and translated to z = 0.2.
  AppendTriangleFan(
      {
          Point3f(0.0f, 0.5f, 0.2f),  //
          Point3f(0.0f, 1.0f, 0.2f),  //
          Point3f(1.0f, 1.0f, 0.2f),  //
          Point3f(1.0f, 0.5f, 0.2f)   //
      },
      {}, &mesh);
  return mesh;
}

TEST(RaytracerTest, TestFindFirstHit) {
  std::unique_ptr<Raytracer> raytracer = Raytracer::Build(MakeTestScene());

  float t_hit;
  int triangle_index;
  bool hit;

  hit = raytracer->FindFirstHit(Point3f(0.30f, 0.25f, -10.0f),
                                Vector3f(0.0f, 0.0f, 2.0f), &t_hit,
                                &triangle_index);
  // The first hit should be at the unit-square in the XY plane, at t = 5, since
  // the ray has z=2.0 and is 10.0 units away from the plane.
  EXPECT_TRUE(hit);
  EXPECT_NEAR(5.0f, t_hit, 1e-5f);

  hit = raytracer->FindFirstHit(Point3f(0.80f, 0.75f, 0.3f),
                                Vector3f(0.0f, 0.0f, -1.0f), &t_hit,
                                &triangle_index);
  // The first hit should be at the rectangle at z=0.2f, at t = 0.1, since
  // the ray has z=-1.0 and is 0.1 units away from the plane.
  EXPECT_TRUE(hit);
  EXPECT_NEAR(0.1, t_hit, 1e-5f);
}

TEST(RaytracerTest, TestCountIntersections) {
  std::unique_ptr<Raytracer> raytracer = Raytracer::Build(MakeTestScene());

  // Note that the query rays specifically do not fall on the diagonal of the
  // quad, since that may result in the raytracer detecting two triangle
  // intersections.
  //
  // Detecting such cases would be expensive, and is unnecessary.

  // No intersections due to t_max.
  EXPECT_EQ(0, raytracer->CountIntersections({0.3f, 0.25f, -2.0f},
                                             {0.0f, 0.0f, 10.0f}, 0.01f, 100));

  EXPECT_EQ(1, raytracer->CountIntersections({0.3f, 0.25f, -2.0f},
                                             {0.0f, 0.0f, 10.0f}, 1.0f, 100));

  EXPECT_EQ(2, raytracer->CountIntersections({0.75f, 0.30f, -2.0f},
                                             {0.0f, 0.0f, 10.0f}, 1.0f, 100));

  EXPECT_EQ(3, raytracer->CountIntersections({0.75f, 0.80f, -2.0f},
                                             {0.0f, 0.0f, 10.0f}, 10.0f, 3));

  // Limited to 1 intersection by max_intersection_count.
  EXPECT_EQ(1, raytracer->CountIntersections({0.75f, 0.81f, -2.0f},
                                             {0.0f, 0.0f, 10.0f}, 1.0f, 1));

  EXPECT_EQ(3, raytracer->CountIntersections(
                   {0.75f, 0.80f, -2.0f}, {0.0f, 0.0f, 10.0f},
                   std::numeric_limits<float>::infinity(),
                   std::numeric_limits<int>::max()));
  EXPECT_EQ(3, raytracer->CountIntersections(
                   {0.75f, 0.80f, -2.0f}, {0.0f, 0.0f, 10.0f},
                   std::numeric_limits<float>::infinity(), 3));
}

TEST(RaytracerTest, TestFindAllIntersections) {
  std::unique_ptr<Raytracer> raytracer = Raytracer::Build(MakeTestScene());

  // No intersections, since the query ray originates far to the left of the
  // scene, with x = -10.
  {
    std::vector<Raytracer::Intersection> intersections;
    raytracer->FindAllIntersections({-10.0f, 0.25f, -2.0f}, {0.0f, 0.0f, 10.0f},
                                    &intersections);
    EXPECT_TRUE(intersections.empty());
  }

  {
    std::vector<Raytracer::Intersection> intersections;
    raytracer->FindAllIntersections({0.3f, 0.25f, -2.0f}, {0.0f, 0.0f, 10.0f},
                                    &intersections);
    EXPECT_EQ(1, intersections.size());
    EXPECT_NEAR(0.2f, intersections[0].t_hit, 1e-5f);
    // This is the second triangle.
    EXPECT_EQ(1, intersections[0].triangle_index);
  }

  {
    std::vector<Raytracer::Intersection> intersections;
    raytracer->FindAllIntersections({0.8f, 0.75f, -2.0f}, {0.0f, 0.0f, 0.1f},
                                    &intersections);
    EXPECT_EQ(3, intersections.size());
  }
}

}  // namespace
}  // namespace geometry
}  // namespace seurat
