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

#ifndef VR_SEURAT_GEOMETRY_QUAD_H_
#define VR_SEURAT_GEOMETRY_QUAD_H_

#include <array>

#include "ion/math/matrix.h"
#include "ion/math/range.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"

namespace seurat {
namespace geometry {

using ion::math::Point2f;

// 2D quad. The vertices must be in counter-clockwise order.
using Quad2i = std::array<ion::math::Point2i, 4>;
using Quad2f = std::array<ion::math::Point2f, 4>;

// 3D quad.
using Quad3f = std::array<ion::math::Point3f, 4>;

// Projects all vertices of the |quad| using the |matrix| and returns them in a
// new quad.
inline Quad3f ProjectQuad(const ion::math::Matrix4f& matrix,
                          const Quad3f& quad) {
  return {{ion::math::ProjectPoint(matrix, quad[0]),
           ion::math::ProjectPoint(matrix, quad[1]),
           ion::math::ProjectPoint(matrix, quad[2]),
           ion::math::ProjectPoint(matrix, quad[3])}};
}

// Flattens a 3D |quad| into a 2D quad by discarding the z-component and returns
// the result.
inline Quad2f FlattenQuad(const Quad3f& quad) {
  return {{Point2f(quad[0][0], quad[0][1]), Point2f(quad[1][0], quad[1][1]),
           Point2f(quad[2][0], quad[2][1]), Point2f(quad[3][0], quad[3][1])}};
}

// Computes the bounding box of a |quad|.
inline ion::math::Range2f QuadBoundingBox(const Quad2f& quad) {
  ion::math::Range2f range;
  for (const Point2f& p : quad) range.ExtendByPoint(p);
  return range;
}

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_QUAD_H_
