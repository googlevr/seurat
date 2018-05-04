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

#ifndef VR_SEURAT_GEOMETRY_TRIANGLE_H_
#define VR_SEURAT_GEOMETRY_TRIANGLE_H_

#include <array>

#include "ion/math/range.h"
#include "ion/math/rangeutils.h"
#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "seurat/geometry/plane.h"

namespace seurat {
namespace geometry {

template <int Dimension, typename T>
using Triangle = std::array<ion::math::Point<Dimension, T>, 3>;

// Returns the unit normal vector of the |triangle|, pointing in the direction
// of the visible halfspace for triangles specified in counter-clockwise order.
template <typename T>
ion::math::Vector<3, T> NormalFromTriangle(const Triangle<3, T>& triangle) {
  return ion::math::Normalized(
      ion::math::Cross(triangle[1] - triangle[0], triangle[2] - triangle[0]));
}

// Returns the the plane of the |triangle|, with positive signed distance to
// points in the visible halfspace for triangles specified in counter-clockwise
// order.
template <typename T>
Plane<3, T> PlaneFromTriangle(const Triangle<3, T>& triangle) {
  return Plane<3, T>(triangle[0], ion::math::Cross(triangle[1] - triangle[0],
                                                   triangle[2] - triangle[0]));
}

// Returns the bounding box of the given |triangle|.
template <int Dimension, typename T>
ion::math::Range<Dimension, T> BoundingBoxFromTriangle(
    const Triangle<Dimension, T>& triangle) {
  ion::math::Range<Dimension, T> bounds;
  for (const auto& v : triangle) {
    bounds.ExtendByPoint(v);
  }
  return bounds;
}

// Returns the barycentric coordinates of the specified point relative to the
// triangle.
template <int Dimension, typename T>
ion::math::Point<3, T> BarycentricFromPoint(
    const Triangle<Dimension, T>& tri,
    const ion::math::Point<Dimension, T>& p) {
  ion::math::Vector<Dimension, T> u = tri[1] - tri[0];
  ion::math::Vector<Dimension, T> v = tri[2] - tri[0];
  ion::math::Vector<Dimension, T> p_tri = p - tri[0];
  T d00 = ion::math::Dot(u, u);
  T d01 = ion::math::Dot(u, v);
  T d11 = ion::math::Dot(v, v);
  T d20 = ion::math::Dot(p_tri, u);
  T d21 = ion::math::Dot(p_tri, v);

  ion::math::Point<3, T> b;
  b[1] = (d11 * d20 - d01 * d21) / (d00 * d11 - d01 * d01);
  b[2] = (d00 * d21 - d01 * d20) / (d00 * d11 - d01 * d01);
  b[0] = T{1} - b[1] - b[2];
  return b;
}

// Converts from barycentric coordinates to a point relative to the triangle.
template <int Dimension, typename T>
ion::math::Point<Dimension, T> PointFromBarycentric(
    const Triangle<Dimension, T>& tri,
    const ion::math::Point<3, T>& barycentric_coords) {
  return tri[0] * barycentric_coords[0] +  //
         tri[1] * barycentric_coords[1] +  //
         tri[2] * barycentric_coords[2];
}

// Returns whether the given barycentric coordinates represent a point within
// the associated triangle.
template <typename T>
bool BarycentricCoordinatesInTriangle(
    const ion::math::Point<3, T>& barycentric_coords) {
  for (int d = 0; d < 3; ++d) {
    if (T(0) > barycentric_coords[d] || barycentric_coords[d] > T(1)) {
      return false;
    }
  }
  return true;
}

using Triangle2f = Triangle<2, float>;
using Triangle3f = Triangle<3, float>;
using Triangle2d = Triangle<2, double>;
using Triangle3d = Triangle<3, double>;

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_TRIANGLE_H_
