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

#include "seurat/tiler/rail_penalty_cost_function.h"

#include <array>
#include <cmath>
#include <vector>

#include "ceres/ceres.h"
#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"

namespace seurat {
namespace tiler {

using ion::math::Vector3d;

bool RailPenaltyCostFunction::Evaluate(double const* const* parameters,
                                       double* residuals,
                                       double** jacobians) const {
  // The plane's normal.
  Vector3d n(parameters[0][0], parameters[0][1], parameters[0][2]);
  ion::math::Normalize(&n);
  // The plane's center.
  Vector3d c(parameters[0][3], parameters[0][4], parameters[0][5]);
  double c_dot_n = ion::math::Dot(c, n);
  const double min_depth = ray_depth_.GetMinPoint();
  const double max_depth = ray_depth_.GetMaxPoint();
  int residual_index = 0;
  for (int corner = 0; corner < 4; ++corner) {
    Vector3d r = Vector3d(rails_[corner]);
    double r_dot_n = ion::math::Dot(r, n);
    double normalized_ray_depth = c_dot_n / r_dot_n;

    residuals[residual_index++] =
        std::max(min_depth - normalized_ray_depth, 0.0);

    residuals[residual_index++] =
        std::max(normalized_ray_depth - max_depth, 0.0);
  }
  for (int i = 0; i < residual_index; ++i) {
    if (!std::isfinite(residuals[i])) {
      return false;
    }
  }
  if (jacobians != nullptr) {
    int jacobian_index = 0;
    for (int corner = 0; corner < 4; ++corner) {
      Vector3d r = Vector3d(rails_[corner]);
      double r_dot_n = ion::math::Dot(r, n);
      double normalized_ray_depth = c_dot_n / r_dot_n;

      // Gradient of corner_min_penalty.
      {
        Vector3d d_corner_min_penalty_d_n = Vector3d::Zero();
        Vector3d d_corner_min_penalty_d_c = Vector3d::Zero();
        if (normalized_ray_depth < min_depth) {
          // Gradient of:
          //   Residual(n) = min_depth - c.n / r.n
          // is
          //   r * c.n/(r.n)^2 - c/r.n
          d_corner_min_penalty_d_n =
              (c_dot_n / (r_dot_n * r_dot_n)) * r - c / r_dot_n;
          // Gradient of:
          //   Residual(c) = min_depth - c.n / r.n
          // is
          //    -n / (r . n)
          d_corner_min_penalty_d_c = -n / r_dot_n;
        }

        jacobians[0][jacobian_index++] = d_corner_min_penalty_d_n[0];
        jacobians[0][jacobian_index++] = d_corner_min_penalty_d_n[1];
        jacobians[0][jacobian_index++] = d_corner_min_penalty_d_n[2];
        jacobians[0][jacobian_index++] = d_corner_min_penalty_d_c[0];
        jacobians[0][jacobian_index++] = d_corner_min_penalty_d_c[1];
        jacobians[0][jacobian_index++] = d_corner_min_penalty_d_c[2];
      }

      // Gradient of corner_max_penalty.
      {
        Vector3d d_corner_max_penalty_d_n = Vector3d::Zero();
        Vector3d d_corner_max_penalty_d_c = Vector3d::Zero();
        if (normalized_ray_depth > max_depth) {
          // Gradient of:
          //   Residual(n) = c.n / r.n - max_depth
          // is
          //   -r * c.n/(r.n)^2 + c/r.n
          d_corner_max_penalty_d_n =
              -(c_dot_n / (r_dot_n * r_dot_n)) * r + c / r_dot_n;
          // Gradient of:
          //   Residual(c) = c.n / r.n - max_depth
          // is
          //   n / (r . n)
          d_corner_max_penalty_d_c = n / r_dot_n;
        }

        jacobians[0][jacobian_index++] = d_corner_max_penalty_d_n[0];
        jacobians[0][jacobian_index++] = d_corner_max_penalty_d_n[1];
        jacobians[0][jacobian_index++] = d_corner_max_penalty_d_n[2];
        jacobians[0][jacobian_index++] = d_corner_max_penalty_d_c[0];
        jacobians[0][jacobian_index++] = d_corner_max_penalty_d_c[1];
        jacobians[0][jacobian_index++] = d_corner_max_penalty_d_c[2];
      }

      // Gradient of corner_max_penalty.
    }
    for (int i = 0; i < jacobian_index; ++i) {
      if (!std::isfinite(jacobians[0][i])) {
        return false;
      }
    }
  }
  return true;
}

}  // namespace tiler
}  // namespace seurat
