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

#ifndef VR_SEURAT_GEOMETRY_POLYGON_H_
#define VR_SEURAT_GEOMETRY_POLYGON_H_

#include <vector>

#include "ion/math/angle.h"
#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/geometry/quad.h"
#include "seurat/image/image.h"

namespace seurat {
namespace geometry {

// Functions related to 2d polygons.

// Extent (range) of polygon in directions |direction| relative to point
// |point|. |direction| is a pair of orthonormal vectors.
ion::math::Range2f Extent(absl::Span<const ion::math::Point2f> polygon,
                          const std::array<ion::math::Vector2f, 2>& direction,
                          const ion::math::Point2f& point);

// Outputs the signed area of |polygon|. The polygon is not required to be
// convex. Here we assume that the polygon is free from self intersections.
// Then the polygon has a positive area when its interior stays on the left side
// while we traverse the sequence of points that form the boundary.
float PolygonArea(absl::Span<const ion::math::Point2f> polygon);

// Outputs the minimum area oriented bounding box that encloses the convex
// polygon |convex_polygon|. The area of the bounding box is also returned.
Quad2f ComputeConvexPolygonOBB(
    absl::Span<const ion::math::Point2f> convex_polygon, float* obb_area);

// Outputs a simple polygon that bounds tightly the non-empty texels in
// |texture|. The output vertices are ordered counter-clockwise.
std::vector<ion::math::Point2f> ComputeTextureBoundary(
    const image::Image4f& texture, int* opaque_texel_count);

// Outputs the angle, in radians, of a counter-clockwise rotation that will make
// the oriented bounding box vertical, with a portrait orientation.
// |pivot_index| is the corner of the oriented bounding box that will become the
// bottom-left corner of the axis-aligned bounding box after rotation.
float ComputeOBBTilt(absl::Span<const ion::math::Point2f> obb,
                     int* pivot_index);

// Returns whether the sequence of points represent a convex counter-clockwise
// loop.
bool IsConvexCounterClockwise(absl::Span<const ion::math::Point2f> points);

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_POLYGON_H_
