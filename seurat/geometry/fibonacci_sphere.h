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

#ifndef VR_SEURAT_GEOMETRY_FIBONACCI_SPHERE_H_
#define VR_SEURAT_GEOMETRY_FIBONACCI_SPHERE_H_

#include "ion/math/vector.h"

namespace seurat {
namespace geometry {

// Generates the i'th point (starting with i=0) on the unit-sphere, from a set
// of |num_points| total points, based on the method described in "Spherical
// Fibonacci Mapping" (Keinert et. al, 2015)
//
// The resulting points are relatively-quick to compute and are close to
// uniformly-sampling the sphere.
//
// |scrambler| should be a random angle uniformly distributed over [0, 2*Pi],
// held constant for all points in a single point set.
ion::math::Point3d GenerateFibonacciSpherePoint(int num_points,
                                                double scrambler, int i);

// Returns the index of the fibonacci-sphere point which is closest to the given
// normalized_direction vector.
int InverseFibonacciSphereMapping(
    int num_points, double scrambler,
    const ion::math::Vector3d& normalized_direction);

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_FIBONACCI_SPHERE_H_
