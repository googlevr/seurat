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

#include "seurat/tiler/plane_projection_cost_function.h"

#include <array>
#include <cmath>
#include <vector>

#include "ceres/ceres.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "seurat/base/util.h"
#include "seurat/geometry/plane.h"

namespace seurat {
namespace tiler {

using ion::math::Vector3d;

bool PlaneProjectionCostFunction::Evaluate(double const* const* parameters,
                                           double* residuals,
                                           double** jacobians) const {
  // The plane's normal.
  Vector3d n(parameters[0][0], parameters[0][1], parameters[0][2]);
  // The plane's center.
  Vector3d c(parameters[0][3], parameters[0][4], parameters[0][5]);
  int residual_index = 0;
  bool has_weights = !point_set_.weights.empty();
  for (int poi : points_of_interest_) {
    Vector3d p(point_set_.positions[poi][0], point_set_.positions[poi][1],
               point_set_.positions[poi][2]);
    double center_dot_normal = ion::math::Dot(c, n);
    double point_dot_normal = ion::math::Dot(p, n);
    // TODO(puneetl):  Optimize away the sqrt by precomputing & caching.
    float weight = has_weights ? std::sqrt(point_set_.weights[poi]) : 1.0f;
    residuals[residual_index] =
        (center_dot_normal / point_dot_normal - 1.0) * weight;
    if (!std::isfinite(residuals[residual_index])) {
      return false;
    }
    residual_index++;
  }
  if (jacobians != nullptr) {
    int jacobian_index = 0;
    for (int poi : points_of_interest_) {
      Vector3d p(point_set_.positions[poi][0], point_set_.positions[poi][1],
                 point_set_.positions[poi][2]);

      const double nx = n[0];
      const double ny = n[1];
      const double nz = n[2];

      const double cx = c[0];
      const double cy = c[1];
      const double cz = c[2];

      const double px = p[0];
      const double py = p[1];
      const double pz = p[2];

      const double p_dot_n = ion::math::Dot(p, n);
      const double p_dot_n_2 = p_dot_n * p_dot_n;

      // The following code was generated via mathematica.

      // d [residual] / d [nx]
      jacobians[0][jacobian_index + 0] =
          (-(cy * ny * px) - cz * nz * px + cx * ny * py + cx * nz * pz) /
          p_dot_n_2;
      // d [residual] / d [ny]
      jacobians[0][jacobian_index + 1] =
          (-((cx * nx + cz * nz) * py) + cy * (nx * px + nz * pz)) / p_dot_n_2;
      // d [residual] / d [nz]
      jacobians[0][jacobian_index + 2] =
          (cz * (nx * px + ny * py) - (cx * nx + cy * ny) * pz) / p_dot_n_2;

      // d [residual] / d [cx]
      jacobians[0][jacobian_index + 3] = nx / p_dot_n;
      // d [residual] / d [cy]
      jacobians[0][jacobian_index + 4] = ny / p_dot_n;
      // d [residual] / d [cz]
      jacobians[0][jacobian_index + 5] = nz / p_dot_n;

      if (!point_set_.weights.empty()) {
        for (int i = 0; i < 6; ++i) {
          // TODO(puneetl):  Optimize this!
          jacobians[0][jacobian_index + i] *=
              std::sqrt(point_set_.weights[poi]);
        }
      }
      for (int i = 0; i < 6; ++i) {
        if (!std::isfinite(jacobians[0][jacobian_index + i])) {
          return false;
        }
      }
      jacobian_index += 6;
    }
  }
  return true;
}

}  // namespace tiler
}  // namespace seurat
