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

#ifndef VR_SEURAT_GEOMETRY_PLANE_H_
#define VR_SEURAT_GEOMETRY_PLANE_H_

#include <algorithm>
#include <array>

#include "ion/math/matrix.h"
#include "ion/math/matrixutils.h"
#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"

namespace seurat {
namespace geometry {

// Stores a plane based on the plane equation x . n + d == 0 as a unit-length
// normal vector, n, and the scalar value, d.
template <int Dimension, typename T>
class Plane {
 public:
  typedef ion::math::Vector<Dimension, T> VectorT;
  typedef ion::math::Point<Dimension, T> PointT;
  typedef Plane<Dimension, T> PlaneT;

  // Default constructor creates an invalid plane, with a zero vector normal.
  Plane() : normal_(VectorT::Zero()), d_(0) {}

  // Create from normal and signed distance to the origin.
  //
  // The normal vector will be normalized internally.
  Plane(const VectorT& normal, T d);

  // Constructs a plane from plane equation coefficients:
  //   a * x + b * y + c * z + d == 0
  //
  // Coefficients will be normalized internally.
  //
  // Only supported when Dimension == 3.
  Plane(T a, T b, T c, T d) : PlaneT({a, b, c}, d) {}

  // Constructs a plane from plane equation coefficients.
  //
  // Coefficients will be normalized internally.
  static PlaneT FromCoefficients(
      const ion::math::Vector<Dimension + 1, T>& coefficients);

  // Create from point on plane and normal vector
  Plane(const PointT& point, const VectorT& normal);

  // Create a plane from the given unit normal vector and signed distance to the
  // origin.
  //
  // No checking or normalization is performed on these values.
  static PlaneT FromNormalizedCoefficients(const VectorT& unit_normal, T d);

  // Returns the unit normal vector, (a, b, c).
  const VectorT& GetNormal() const { return normal_; }

  // Returns the signed distance from the origin to the plane in the direction
  // opposite the normal.
  T GetD() const { return d_; }

  // Returns all plane-equation coefficients as a single vector.
  ion::math::Vector<Dimension + 1, T> GetCoefficients() const;

  // Returns true if all values in two instances are equal.
  static bool AreValuesEqual(const Plane& p0, const Plane& p1);

  // Returns false if the plane is degenerate.
  bool IsValid() const;

  // Returns the signed distance from the plane to the point.
  //
  // Positive values are in the direction of the normal.
  T SignedDistanceToPoint(const PointT& point) const;

  // Returns the projection of the |point| on the plane, along the normal
  // vector.
  PointT ProjectPoint(const PointT& point) const;

  // Returns the 4x4 projection matrix that projects 3D points onto this plane
  // using the given |center_of_projection|. If the |center_of_projection| is in
  // the negative half-space of the plane, projected points will have a negative
  // w-component, if they are 'behind' the center of projection. 'behind' means
  // that the ray from the center of projection through the point (before
  // projection) has no valid (positive t-value) intersection with the plane.
  //
  // This is only implemented for Dimension = 3.
  ion::math::Matrix<Dimension + 1, T> ProjectionMatrix(
      const PointT& center_of_projection) const;

  // Intersects the plane with a ray defined by |origin| and |direction|. If a
  // valid intersection is found (t > 0), the method returns true and stores the
  // t value at the intersection point in |t_hit|. Otherwise the method returns
  // false and does not modify |t_hit|. Rays with the origin in the plane don't
  // intersect it.
  bool IntersectRay(const PointT& origin, const VectorT& direction,
                    T* t_hit) const;

  // Returns the plane with reverse orientation.  The plane's geometry is the
  // same, but the normal orientation is in the opposite direction.
  PlaneT GetReversePlane() const;

  // Returns the plane transformed with the given matrix.
  //
  // Note:  This does NOT compute the inverse-transpose of the given |matrix|,
  // which is necessary to transform a plane correctly.
  PlaneT Transform(
      const ion::math::Matrix<Dimension + 1, T>& normal_matrix) const;

  // Returns a normalized vector tangent to the plane.
  VectorT GetTangent() const;

 private:
  VectorT normal_;
  T d_;
};

template <int Dimension, typename T>
Plane<Dimension, T>::Plane(const VectorT& normal, T d) {
  const auto length_squared = ion::math::LengthSquared(normal);
  DCHECK_GT(length_squared, static_cast<T>(0));
  if (length_squared == static_cast<T>(1)) {
    // If the given normal is already unit-length, then use it.
    normal_ = normal;
    d_ = d;
  } else {
    // Normalize
    const auto length = std::sqrt(length_squared);
    normal_ = normal / length;
    d_ = d / length;
  }
}

template <int Dimension, typename T>
Plane<Dimension, T>::Plane(const PointT& point, const VectorT& normal) {
  normal_ = ion::math::Normalized(normal);
  d_ = -ion::math::Dot(normal_, point - PointT::Zero());
}

template <int Dimension, typename T>
Plane<Dimension, T> Plane<Dimension, T>::FromNormalizedCoefficients(
    const VectorT& unit_normal, T d) {
  PlaneT plane;
  plane.normal_ = unit_normal;
  plane.d_ = d;
  return plane;
}

template <int Dimension, typename T>
Plane<Dimension, T> Plane<Dimension, T>::FromCoefficients(
    const ion::math::Vector<Dimension + 1, T>& coefficients) {
  VectorT normal;
  for (int d = 0; d < Dimension; ++d) {
    normal[d] = coefficients[d];
  }
  return PlaneT(normal, coefficients[Dimension]);
}

template <int Dimension, typename T>
ion::math::Vector<Dimension + 1, T> Plane<Dimension, T>::GetCoefficients()
    const {
  return ion::math::Vector<Dimension + 1, T>(GetNormal(), GetD());
}

template <int Dimension, typename T>
bool Plane<Dimension, T>::AreValuesEqual(const PlaneT& p0, const PlaneT& p1) {
  return VectorT::AreValuesEqual(p0.normal_, p1.normal_) && p0.d_ == p1.d_;
}

template <int Dimension, typename T>
T Plane<Dimension, T>::SignedDistanceToPoint(const PointT& point) const {
  return ion::math::Dot(point - PointT::Zero(), normal_) + d_;
}

template <int Dimension, typename T>
bool Plane<Dimension, T>::IsValid() const {
  return normal_ != VectorT::Zero();
}

template <int Dimension, typename T>
ion::math::Point<Dimension, T> Plane<Dimension, T>::ProjectPoint(
    const PointT& point) const {
  const T distance = SignedDistanceToPoint(point);
  return point - normal_ * distance;
}

template <int Dimension, typename T>
ion::math::Matrix<Dimension + 1, T> Plane<Dimension, T>::ProjectionMatrix(
    const PointT& center_of_projection) const {
  // Flip the plane's orientation, if the |center_of_projection| is in the
  // positive half-space of the plane. Otherwise, points 'behind' the
  // |center_of_projection| will have a positive w-component after the
  // projection.
  if (SignedDistanceToPoint(center_of_projection) > static_cast<T>(0)) {
    return GetReversePlane().ProjectionMatrix(center_of_projection);
  }

  const PointT& c = center_of_projection;
  const VectorT& n = normal_;
  const T dist = SignedDistanceToPoint(c);
  // For a derivation of the matrix see:
  // "Matrix of projection on a plane"
  // http://maverick.inria.fr/~Xavier.Decoret/resources/maths/plane-projection.pdf
  // clang-format off
  return ion::math::Matrix<4, T>(
      n[0]*c[0] - dist, n[1]*c[0],        n[2]*c[0],        d_*c[0],
      n[0]*c[1],        n[1]*c[1] - dist, n[2]*c[1],        d_*c[1],
      n[0]*c[2],        n[1]*c[2],        n[2]*c[2] - dist, d_*c[2],
      n[0],             n[1],             n[2],             d_ - dist);
  // clang-format on
}

template <int Dimension, typename T>
bool Plane<Dimension, T>::IntersectRay(const PointT& origin,
                                       const VectorT& direction,
                                       T* t_hit) const {
  // Check if the ray is parallel to the plane.
  const T normal_dot_direction = ion::math::Dot(normal_, direction);
  if (normal_dot_direction == static_cast<T>(0)) return false;

  // Compute the t value at the intersection point and store it in |t_hit| if it
  // is greater than zero.
  const T t = -SignedDistanceToPoint(origin) / normal_dot_direction;
  const bool hit = (t > static_cast<T>(0));
  if (hit) *t_hit = t;

  return hit;
}

template <int Dimension, typename T>
Plane<Dimension, T> Plane<Dimension, T>::GetReversePlane() const {
  return PlaneT::FromNormalizedCoefficients(-normal_, -d_);
}

template <int Dimension, typename T>
Plane<Dimension, T> Plane<Dimension, T>::Transform(
    const ion::math::Matrix<Dimension + 1, T>& normal_matrix) const {
  const ion::math::Vector<Dimension + 1, T> plane_coeffs(normal_, d_);
  const auto transformed_coeffs = normal_matrix * plane_coeffs;
  VectorT normal;
  for (int d = 0; d < Dimension; ++d) {
    normal[d] = transformed_coeffs[d];
  }
  return PlaneT(normal, transformed_coeffs[Dimension]);
}

template <int Dimension, typename T>
ion::math::Vector<Dimension, T> Plane<Dimension, T>::GetTangent() const {
  // Select the unit basis vector for which the projection onto normal_ has the
  // smallest magnitude.
  int smallest_normal_component_index = 0;
  T smallest_normal_component_value = ion::math::Abs(normal_[0]);
  for (int d = 1; d < Dimension; ++d) {
    if (ion::math::Abs(normal_[d]) < smallest_normal_component_value) {
      smallest_normal_component_value = ion::math::Abs(normal_[d]);
      smallest_normal_component_index = d;
    }
  }

  VectorT basis_vector = VectorT::Zero();
  basis_vector[smallest_normal_component_index] = T{1};

  // Orthogonalize this vector with respect to normal_.
  //
  // Instead of naively computing n * (n . b) / (n . n), we can simplify:
  //  * n . n == 1 since normal_ is stored as a unit-vector.
  //  * n . b is simply the element of n with index corresponding to the
  //  nonzero element of the basis vector, b.
  VectorT projection = normal_ * normal_[smallest_normal_component_index];
  return ion::math::Normalized(basis_vector - projection);
}

typedef Plane<2, float> Plane2f;
typedef Plane<2, double> Plane2d;
typedef Plane<3, float> Plane3f;
typedef Plane<3, double> Plane3d;

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_PLANE_H_
