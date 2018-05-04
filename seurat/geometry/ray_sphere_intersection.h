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

#ifndef VR_SEURAT_GEOMETRY_RAY_SPHERE_INTERSECTION_H_
#define VR_SEURAT_GEOMETRY_RAY_SPHERE_INTERSECTION_H_

#include "ion/math/vector.h"

namespace seurat {
namespace geometry {

// Finds the first ray-sphere intersection, returning the t-value associated
// with the intersection point.
//
// Returns whether there is an intersection.
bool ComputeRaySphereIntersection(const ion::math::Point3f& sphere_center,
                                  float sphere_radius,
                                  const ion::math::Point3f& ray_start,
                                  const ion::math::Vector3f& ray_direction,
                                  float* t_hit);

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_RAY_SPHERE_INTERSECTION_H_
