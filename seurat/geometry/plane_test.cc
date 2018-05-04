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

#include "seurat/geometry/plane.h"

#include "ion/math/matrix.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "gtest/gtest.h"

namespace seurat {
namespace geometry {
namespace {

using ion::math::Matrix4f;
using ion::math::Point3d;
using ion::math::Point3f;
using ion::math::Vector3d;
using ion::math::Vector3f;
using ion::math::Vector4f;

TEST(PlaneTest, ValidInvalid) {
  const Point3f point(0.0f, 0.0f, 0.0f);
  const Vector3f normal(0.0f, 0.0f, 1.0f);

  const Plane3f valid_plane(point, normal);
  EXPECT_TRUE(valid_plane.IsValid());

  const Plane3f invalid_plane;
  EXPECT_FALSE(invalid_plane.IsValid());
}

TEST(PlaneTest, CreateFromPointAndNormal) {
  const Plane3f plane =
      Plane3f(Point3f(3.0f, 1.0f, 1.0f), Vector3f(1.0f, 0.0f, 0.0f));
  const auto& normal = plane.GetNormal();

  EXPECT_FLOAT_EQ(1.0f, ion::math::LengthSquared(normal));
  EXPECT_FLOAT_EQ(1.0f, normal[0]);
  EXPECT_FLOAT_EQ(0.0f, normal[1]);
  EXPECT_FLOAT_EQ(0.0f, normal[2]);
  EXPECT_FLOAT_EQ(-3.0f, plane.GetD());
}

TEST(PlaneTest, FromNormalizedCoefficients) {
  // Test with a well-behaved function call
  const auto unit_normal = ion::math::Normalized(Vector3f(2.0f, -2.0f, 4.0f));
  const auto plane = Plane3f::FromNormalizedCoefficients(unit_normal, -6.0f);

  EXPECT_FLOAT_EQ(unit_normal[0], plane.GetNormal()[0]);
  EXPECT_FLOAT_EQ(unit_normal[1], plane.GetNormal()[1]);
  EXPECT_FLOAT_EQ(unit_normal[2], plane.GetNormal()[2]);
  EXPECT_FLOAT_EQ(-6.0f, plane.GetD());

  // No normalization should be performed (for the sake of speed) even if the
  // input is incorrectly unnormalized.
  const auto unnormalized_plane =
      Plane3f::FromNormalizedCoefficients(Vector3f(2.0f, -2.0f, 4.0f), -6.0f);
  EXPECT_FLOAT_EQ(2.0f, unnormalized_plane.GetNormal()[0]);
  EXPECT_FLOAT_EQ(-2.0f, unnormalized_plane.GetNormal()[1]);
  EXPECT_FLOAT_EQ(4.0f, unnormalized_plane.GetNormal()[2]);
  EXPECT_FLOAT_EQ(-6.0f, unnormalized_plane.GetD());
}

TEST(PlaneTest, CreateFromUnnormalizedCoefficients) {
  const Plane3f plane{5.0f, 0.0f, 0.0f, 10.0f};

  EXPECT_FLOAT_EQ(1.0f, ion::math::LengthSquared(plane.GetNormal()));
  EXPECT_FLOAT_EQ(2.0f,
                  plane.SignedDistanceToPoint(Point3f(0.0f, -42.0f, 37.0f)));
  EXPECT_FLOAT_EQ(2.0f, plane.GetD());
  EXPECT_TRUE(Plane3f::AreValuesEqual(
      plane,
      Plane3f(Point3f(-2.0f, -99.0f, 99.0f), Vector3f(100.0f, 0.0f, 0.0f))));
}

TEST(PlaneTest, CreateFromUnnormalizedCoefficientsVector) {
  const Plane3f plane =
      Plane3f::FromCoefficients(Vector4f(5.0f, 0.0f, 0.0f, 10.0f));

  EXPECT_FLOAT_EQ(1.0f, ion::math::LengthSquared(plane.GetNormal()));
  EXPECT_FLOAT_EQ(2.0f,
                  plane.SignedDistanceToPoint(Point3f(0.0f, -42.0f, 37.0f)));
  EXPECT_FLOAT_EQ(2.0f, plane.GetD());
  EXPECT_TRUE(Plane3f::AreValuesEqual(
      plane,
      Plane3f(Point3f(-2.0f, -99.0f, 99.0f), Vector3f(100.0f, 0.0f, 0.0f))));

  EXPECT_EQ(Vector4f(1.0f, 0.0f, 0.0f, 2.0f), plane.GetCoefficients());
}

TEST(PlaneTest, CreateFromPointAndNormalUnnormalized) {
  // Create a plane with a point and an unnormalized vector
  const Point3f point(3.0f, 1.0f, 1.0f);
  const Vector3f direction(3.0f, 3.0f, 0.0f);

  const Plane3f plane = Plane3f(point, direction);

  // The normal should be unit-length.
  EXPECT_FLOAT_EQ(1.0f, ion::math::LengthSquared(plane.GetNormal()));
  // The normal should be in the correct direction.
  EXPECT_FLOAT_EQ(1.0f, ion::math::Dot(ion::math::Normalized(direction),
                                       plane.GetNormal()));
  // The point should be on the plane.
  EXPECT_FLOAT_EQ(0.0f, plane.SignedDistanceToPoint(point));
}

TEST(PlaneTest, AreValuesEqual) {
  const Plane3f plane(Point3f(3.0f, 0.0f, 0.0f), Vector3f(1.0f, 0.0f, 0.0f));

  EXPECT_TRUE(Plane3f::AreValuesEqual(plane, Plane3f(1.0f, 0.0f, 0.0f, -3.0f)));

  EXPECT_FALSE(
      Plane3f::AreValuesEqual(plane, Plane3f(1.0f, 0.01f, 0.0f, -3.0f)));
  EXPECT_FALSE(
      Plane3f::AreValuesEqual(plane, Plane3f(1.0f, 0.0f, 0.0f, -3.01f)));
}

TEST(PlaneTest, Sidedness) {
  const Plane3f plane(Point3f(1.0f, 1.0f, 1.0f), Vector3f(1.0f, 1.0f, 1.0f));

  EXPECT_LT(0.0f, plane.SignedDistanceToPoint(Point3f(2.0f, 2.0f, 2.0f)));
  EXPECT_GT(0.0f, plane.SignedDistanceToPoint(Point3f::Zero()));

  const Plane3f reverse = plane.GetReversePlane();

  EXPECT_GT(0.0f, reverse.SignedDistanceToPoint(Point3f(2.0f, 2.0f, 2.0f)));
  EXPECT_LT(0.0f, reverse.SignedDistanceToPoint(Point3f::Zero()));
}

TEST(PlaneTest, Distance) {
  const Plane3f plane(Point3f(3.0f, 0.0f, 0.0f),
                      ion::math::Normalized(Vector3f(1.0f, 0.0f, 0.0f)));

  EXPECT_EQ(0.0f, plane.SignedDistanceToPoint(Point3f(3.0f, 1.0f, 1.0f)));
  EXPECT_EQ(-3.0f, plane.SignedDistanceToPoint(Point3f::Zero()));
}

TEST(PlaneTest, TransformedPlane) {
  const Point3d point(1.0f, 1.0f, 1.0f);
  const Vector3d direction(2.0f, 1.0f, 1.0f);
  const Plane3d plane(point, direction);

  auto transform = ion::math::TranslationMatrix(Point3d(5.0f, 3.0f, 9.0f)) *
                   ion::math::ScaleMatrixH(Vector3d(1.0f, 0.5f, 3.0f));

  const auto transformed_plane =
      plane.Transform(ion::math::Inverse(ion::math::Transpose(transform)));

  EXPECT_DOUBLE_EQ(1.0,
                   ion::math::LengthSquared(transformed_plane.GetNormal()));
  EXPECT_DOUBLE_EQ(0.0,
                   transformed_plane.SignedDistanceToPoint(transform * point));

  const auto point_on_positive_side = point + direction * 3;
  EXPECT_LT(0.0, transformed_plane.SignedDistanceToPoint(
                     transform * point_on_positive_side));
}

TEST(PlaneTest, ProjectPoint) {
  const Plane3f plane(Point3f(3.0f, 0.0f, 0.0f),
                      ion::math::Normalized(Vector3f(1.0f, 0.0f, 0.0f)));
  EXPECT_EQ(Point3f(3.0f, 2.0f, -1.0f),
            plane.ProjectPoint(Point3f(4.3f, 2.0f, -1.0f)));
  EXPECT_EQ(Point3f(3.0f, -1.0f, 4.0f),
            plane.ProjectPoint(Point3f(-1.1f, -1.0f, 4.0f)));
}

TEST(PlaneTest, IntersectRay) {
  const Plane3f plane(Point3f(3.0f, 0.0f, 0.0f),
                      ion::math::Normalized(Vector3f(1.0f, 0.0f, 0.0f)));
  bool hit;
  float t_hit;

  Point3f origin(1.0f, 4.0f, -1.0f);
  Vector3f direction(1.0f, 0.0f, 0.0f);
  hit = plane.IntersectRay(origin, direction, &t_hit);
  EXPECT_TRUE(hit);
  EXPECT_EQ(2.0f, t_hit);

  // Ray pointing away from plane doesn't intersect.
  origin = Point3f(5.0, 2.0f, -1.0f);
  direction = Vector3f(1.0f, 0.0f, 0.0f);
  hit = plane.IntersectRay(origin, direction, &t_hit);
  EXPECT_FALSE(hit);

  // Ray pointing towards plane does intersect.
  origin = Point3f(5.0, 2.0f, -1.0f);
  direction = Vector3f(-1.0f, 0.0f, 0.0f);
  hit = plane.IntersectRay(origin, direction, &t_hit);
  EXPECT_TRUE(hit);
  EXPECT_EQ(2.0f, t_hit);

  // Ray parallel to plane doesn't intersect.
  origin = Point3f(5.0f, 2.0f, -3.0f);
  direction = Vector3f(0.0f, 1.0f, 0.0f);
  hit = plane.IntersectRay(origin, direction, &t_hit);
  EXPECT_FALSE(hit);

  // Ray with origin in plane doesn't intersect.
  origin = Point3f(3.0, 2.0f, -1.0f);
  direction = Vector3f(-1.0f, 0.0f, 0.0f);
  hit = plane.IntersectRay(origin, direction, &t_hit);
  EXPECT_FALSE(hit);
}

TEST(PlaneTest, GetTangentForAxisAlignedPlane) {
  const Plane3f plane(Point3f(3.0f, 0.0f, 0.0f),
                      ion::math::Normalized(Vector3f(1.0f, 0.0f, 0.0f)));

  const auto tangent = plane.GetTangent();
  EXPECT_GT(1e-6f, std::fabs(ion::math::LengthSquared(tangent) - 1.0f));
  EXPECT_GT(1e-6f, ion::math::Dot(tangent, plane.GetNormal()));
}

TEST(PlaneTest, GetTangentForPlaneWithNegativeNormalComponent) {
  const Plane3f plane(Point3f(0.0f, 0.0f, 0.0f),
                      ion::math::Normalized(Vector3f(0.0f, 0.0f, -1.0f)));

  const auto tangent = plane.GetTangent();
  EXPECT_GT(1e-6f, std::fabs(ion::math::LengthSquared(tangent) - 1.0f));
  EXPECT_GT(1e-6f, ion::math::Dot(tangent, plane.GetNormal()));
}

TEST(PlaneTest, GetTangentForNonAxisAlignedPlane) {
  const Plane3d plane(Point3d(3.0, 0.0, 0.0),
                      ion::math::Normalized(Vector3d(2.0, 3.0, 1.6)));

  const auto tangent = plane.GetTangent();
  EXPECT_GT(1e-6f, std::fabs(ion::math::LengthSquared(tangent) - 1.0f));
  EXPECT_GT(1e-6f, ion::math::Dot(tangent, plane.GetNormal()));
}

TEST(PlaneTest, ProjectionMatrix) {
  // The center of projection is in the negative half-space of the plane. The
  // w-component of projected points will be positive, if the ray from the
  // center of projection through the point (before projection) has a valid
  // intersection with the plane.
  const Plane3f plane(Point3f(0.0f, 0.0f, -10.0f),
                      ion::math::Normalized(Vector3f(-0.1f, -0.1f, -1.0f)));
  const Point3f center_of_projection(1.0, 1.0f, 0.0f);
  const Matrix4f projection = plane.ProjectionMatrix(center_of_projection);
  // Test a few points and check if their projections are on the plane.
  std::vector<Point3f> points_world{
      Point3f(1.0f, 1.0f, -1.0f),   // between plane and center of projection
      Point3f(0.0f, 0.0f, -5.0f),   // between plane and center of projection
      Point3f(1.0f, 1.0f, 1.0f),    // behind the center of projection
      Point3f(1.0f, 1.0f, -20.0f),  // on the other side of the plane
      Point3f(0.0f, 0.0f, -10.0f)   // on the plane
  };
  // Double check that the selected points are on the expected side of the
  // plane.
  EXPECT_LT(plane.SignedDistanceToPoint(points_world[0]), 0.0f);
  EXPECT_LT(plane.SignedDistanceToPoint(points_world[1]), 0.0f);
  EXPECT_LT(plane.SignedDistanceToPoint(points_world[2]), 0.0f);
  EXPECT_GT(plane.SignedDistanceToPoint(points_world[3]), 0.0f);
  EXPECT_FLOAT_EQ(plane.SignedDistanceToPoint(points_world[4]), 0.0f);

  for (const Point3f& p_world : points_world) {
    // Compute projected point and check that it is on the plane.
    const Point3f p_plane = ion::math::ProjectPoint(projection, p_world);
    EXPECT_NEAR(0.0f, plane.SignedDistanceToPoint(p_plane), 1e-6f);

    // Compute the projected point in homogenous coordinates.
    const Vector4f p_plane_hom =
        projection * Vector4f(p_world - Point3f::Zero(), 1.0f);

    // Intersect a ray from the center of projection through the point with the
    // plane.
    const Vector3f ray_dir = p_world - center_of_projection;
    float t_hit;
    const bool hit = plane.IntersectRay(center_of_projection, ray_dir, &t_hit);
    if (hit) {
      // If we find a valid intersection, the w-component should be positive and
      // the intersection point should be identical to the projected point.
      EXPECT_GT(p_plane_hom[3], 0.0f);
      const Point3f hit_point = center_of_projection + t_hit * ray_dir;
      ion::math::PointsAlmostEqual(hit_point, p_plane, 1e-6f);
    } else {
      // If the ray misses the plane (i.e. intersects at a negative t-value),
      // the w-component should be negative.
      EXPECT_LT(p_plane_hom[3], 0.0f);
    }
  }
}

}  // namespace
}  // namespace geometry
}  // namespace seurat
