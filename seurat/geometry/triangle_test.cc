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

#include "seurat/geometry/triangle.h"

#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "gtest/gtest.h"
#include "seurat/geometry/plane.h"

namespace seurat {
namespace geometry {
namespace {

using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Vector3f;

TEST(TriangleTest, NormalFromTriangle) {
  // A triangle in the xy plane, in counter-clockwise order when looking in the
  // negative-z direction.
  Triangle3f triangle = {{{0.0f, 0.0f, 0.0f},  //
                          {2.0f, 1.0f, 0.0f},  //
                          {-1.0f, 3.0f, 0.0f}}};

  Vector3f normal = NormalFromTriangle(triangle);

  EXPECT_EQ(Vector3f(0.0f, 0.0f, 1.0f), normal);
  EXPECT_NEAR(1.0f, ion::math::Length(normal), 1.0e-6f);
}

TEST(TriangleTest, PlaneFromTriangle) {
  // A triangle approximately in the xy plane, in counter-clockwise order when
  // looking in the -z direction.
  Triangle3f triangle = {{{0.0f, 0.0f, 0.1f},   //
                          {2.0f, 1.0f, -0.1f},  //
                          {-1.0f, 3.0f, 0.0f}}};

  Plane3f plane = PlaneFromTriangle(triangle);

  EXPECT_NEAR(1.0f, ion::math::Length(plane.GetNormal()), 1.0e-6f);
  EXPECT_LT(0.5f, plane.GetNormal()[2]);
  EXPECT_NEAR(0.0f, plane.SignedDistanceToPoint(triangle[0]), 1.0e-6f);
  EXPECT_NEAR(0.0f, plane.SignedDistanceToPoint(triangle[1]), 1.0e-6f);
  EXPECT_NEAR(0.0f, plane.SignedDistanceToPoint(triangle[2]), 1.0e-6f);

  const Point3f point_in_visible_halfspace = Point3f(0.0f, 0.0f, 10.0f);
  EXPECT_LT(0.0f, plane.SignedDistanceToPoint(point_in_visible_halfspace));
}

TEST(TriangleTest, BarycentricFromPoint) {
  // A triangle almost-parallel to the xy plane, in counter-clockwise order
  // when looking down.
  Triangle3f triangle = {{{0.0f, 0.0f, 0.47f},  //
                          {2.0f, 1.0f, 0.5f},   //
                          {-1.0f, 3.0f, 0.55f}}};

  Point3f point_in_triangle =
      triangle[0] * 0.5f + triangle[1] * 0.3f + triangle[2] * 0.2f;

  Point3f barycentric_coords =
      BarycentricFromPoint(triangle, point_in_triangle);

  EXPECT_NEAR(0.5f, barycentric_coords[0], 1.0e-6f);
  EXPECT_NEAR(0.3f, barycentric_coords[1], 1.0e-6f);
  EXPECT_NEAR(0.2f, barycentric_coords[2], 1.0e-6f);

  EXPECT_EQ(Point3f(1.0f, 0.0f, 0.0f),
            BarycentricFromPoint(triangle, triangle[0]));
  EXPECT_EQ(Point3f(0.0f, 1.0f, 0.0f),
            BarycentricFromPoint(triangle, triangle[1]));
  EXPECT_EQ(Point3f(0.0f, 0.0f, 1.0f),
            BarycentricFromPoint(triangle, triangle[2]));
}

TEST(TriangleTest, BarycentricCoordinatesInTriangle) {
  // A triangle and 3 points, 'a', 'b', 'c', outside the triangle, and point 'd'
  // inside the triangle.
  //
  //   ^
  //   +
  //   |\
  //   | \
  //  a|  \ b
  //   | d \
  //   |    \
  //  -+-----+-->
  //      c

  Triangle2f triangle = {{{0.0f, 0.0f},  //
                          {1.0f, 0.0f},  //
                          {0.0f, 1.0f}}};

  Point2f a(0.51, 0.51);
  Point2f b(0.5, -0.1);
  Point2f c(-0.1, 0.5);
  Point2f d(0.25, 0.25);

  EXPECT_FALSE(
      BarycentricCoordinatesInTriangle(BarycentricFromPoint(triangle, a)));
  EXPECT_FALSE(
      BarycentricCoordinatesInTriangle(BarycentricFromPoint(triangle, b)));
  EXPECT_FALSE(
      BarycentricCoordinatesInTriangle(BarycentricFromPoint(triangle, c)));
  EXPECT_TRUE(
      BarycentricCoordinatesInTriangle(BarycentricFromPoint(triangle, d)));

  // The vertices of the triangle should be within the triangle.
  EXPECT_TRUE(BarycentricCoordinatesInTriangle(
      BarycentricFromPoint(triangle, triangle[0])));
  EXPECT_TRUE(BarycentricCoordinatesInTriangle(
      BarycentricFromPoint(triangle, triangle[1])));
  EXPECT_TRUE(BarycentricCoordinatesInTriangle(
      BarycentricFromPoint(triangle, triangle[2])));
  // A point on the edge should be within the triangle.
  EXPECT_TRUE(BarycentricCoordinatesInTriangle(
      BarycentricFromPoint(triangle, 0.5f * (triangle[1] + triangle[0]))));
}

TEST(TriangleTest, PointFromBarycentric) {
  // A triangle almost-parallel to the xy plane, in counter-clockwise order
  // when looking down.
  Triangle3f triangle = {{{0.0f, 0.0f, 0.47f},  //
                          {2.0f, 1.0f, 0.5f},   //
                          {-1.0f, 3.0f, 0.55f}}};

  Point3f point_in_triangle =
      triangle[0] * 0.5f + triangle[1] * 0.3f + triangle[2] * 0.2f;

  EXPECT_EQ(point_in_triangle,
            PointFromBarycentric(
                triangle, BarycentricFromPoint(triangle, point_in_triangle)));

  EXPECT_EQ(triangle[0],
            PointFromBarycentric(triangle, Point3f(1.0f, 0.0f, 0.0f)));
  EXPECT_EQ(triangle[1],
            PointFromBarycentric(triangle, Point3f(0.0f, 1.0f, 0.0f)));
  EXPECT_EQ(triangle[2],
            PointFromBarycentric(triangle, Point3f(0.0f, 0.0f, 1.0f)));
}

}  // namespace
}  // namespace geometry
}  // namespace seurat
