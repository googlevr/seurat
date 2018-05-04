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

#include "seurat/geometry/ray_sphere_intersection.h"

#include "ion/math/vectorutils.h"

namespace seurat {
namespace geometry {

bool ComputeRaySphereIntersection(const ion::math::Point3f& sphere_center,
                                  float sphere_radius,
                                  const ion::math::Point3f& ray_start,
                                  const ion::math::Vector3f& ray_direction,
                                  float* t_hit) {
  float a = ion::math::LengthSquared(ray_direction);
  float b = 2.0f * ion::math::Dot(ray_start - sphere_center, ray_direction);
  float c = ion::math::LengthSquared(ray_start - sphere_center) -
            (sphere_radius * sphere_radius);
  float discriminant = b * b - 4 * a * c;
  if (discriminant < 0.0f || a == 0.0f) {
    return false;
  }
  float sqrt_discriminant = std::sqrt(discriminant);
  float t0 = (-b + sqrt_discriminant) / (2.0f * a);
  float t1 = (-b - sqrt_discriminant) / (2.0f * a);
  float t_low = std::min(t0, t1);
  float t_high = std::max(t0, t1);
  if (t_low >= 0.0f) {
    *t_hit = t_low;
    return true;
  }
  if (t_high >= 0.0f) {
    *t_hit = t_high;
    return true;
  }
  return false;
}

}  // namespace geometry
}  // namespace seurat
