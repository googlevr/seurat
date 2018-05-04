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

#ifndef VR_SEURAT_GEOMETRY_CONVEX_HULL2D_H_
#define VR_SEURAT_GEOMETRY_CONVEX_HULL2D_H_

#include <vector>
#include "ion/math/vector.h"

namespace seurat {
namespace geometry {

// Returns the convex hull of the given |points|. |points| is in general an
// unorganized set of points rather than the sequence of vertices of some
// polygon. Implements Andrew's monotone chain convex hull algorithm. In the
// output the vertices of the convex hull are enumerated in counter-clockwise
// order.
std::vector<ion::math::Point2f> ComputeConvexHull(
    const std::vector<ion::math::Point2f>& points);

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_CONVEX_HULL2D_H_
