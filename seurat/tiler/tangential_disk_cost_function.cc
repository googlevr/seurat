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

#include "seurat/tiler/tangential_disk_cost_function.h"

#include <array>
#include <cmath>
#include <vector>

#include "ceres/ceres.h"
#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"

namespace seurat {
namespace tiler {

using ion::math::Vector3d;

namespace {

// To be able to directly copy-paste the output of Mathematica's generated
// jacobian expressions directly into code without modification (which is
// convenient if/when the cost function is tweaked), we must define
// Power and Sqrt here.

double Power(double v, int p) {
  // For the sake of efficiency, special-case for the powers we need with the
  // hope that the compiler will inline the multiplication.
  if (p == 2) {
    return v * v;
  } else if (p == 3) {
    return v * v * v;
  } else {
    LOG(FATAL) << "Not supported";
    return 0.0;
  }
}

double Power(double v, double p) {
  // For the sake of efficiency, special-case for the powers we need.
  if (p == 1.5) {
    return std::sqrt(v) * v;
  } else {
    LOG(FATAL) << "Not supported";
    return 0.0;
  }
}

double Sqrt(double d) { return std::sqrt(d); }

}  // namespace

bool TangentialDiskCostFunction::Evaluate(double const* const* parameters,
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
    float weight = has_weights ? std::sqrt(point_set_.weights[poi]) : 1.0f;
    residuals[residual_index] =
        weight *
        ion::math::Length(c - p * (center_dot_normal / point_dot_normal)) /
        ion::math::Length(c);
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

      const double c_norm_2 = ion::math::LengthSquared(c);
      const double c_norm = std::sqrt(c_norm_2);
      const double c_dot_n = ion::math::Dot(c, n);
      const double p_dot_n = ion::math::Dot(p, n);
      const double p_dot_n_2 = p_dot_n * p_dot_n;
      const double p_dot_n_3 = p_dot_n_2 * p_dot_n;

      const double tangential_dist_2 =
          ion::math::LengthSquared(c - p * (c_dot_n / p_dot_n));
      const double tangential_dist = std::sqrt(tangential_dist_2);

      // The following code was generated via mathematica.

      jacobians[0][jacobian_index + 0] =
          ((-(cy * ny * px) - cz * nz * px + cx * ny * py + cx * nz * pz) *
           (cz * (nz * (Power(px, 2) + Power(py, 2)) -
                  (nx * px + ny * py) * pz) +
            cy * (-(py * (nx * px + nz * pz)) +
                  ny * (Power(px, 2) + Power(pz, 2))) +
            cx * (-(ny * px * py) - nz * px * pz +
                  nx * (Power(py, 2) + Power(pz, 2))))) /
          (c_norm * p_dot_n_3 * tangential_dist);
      // d ["tangential cost"] / d [ny]
      jacobians[0][jacobian_index + 1] =
          -((-((cx * nx + cz * nz) * py) + cy * (nx * px + nz * pz)) *
            (cz * (-(nz * (Power(px, 2) + Power(py, 2))) +
                   (nx * px + ny * py) * pz) +
             cy * (py * (nx * px + nz * pz) -
                   ny * (Power(px, 2) + Power(pz, 2))) +
             cx * (ny * px * py + nz * px * pz -
                   nx * (Power(py, 2) + Power(pz, 2))))) /
          (c_norm * p_dot_n_3 * tangential_dist);
      // d ["tangential cost"] / d [nz]
      jacobians[0][jacobian_index + 2] =
          -((cz * (nx * px + ny * py) - (cx * nx + cy * ny) * pz) *
            (cz * (-(nz * (Power(px, 2) + Power(py, 2))) +
                   (nx * px + ny * py) * pz) +
             cy * (py * (nx * px + nz * pz) -
                   ny * (Power(px, 2) + Power(pz, 2))) +
             cx * (ny * px * py + nz * px * pz -
                   nx * (Power(py, 2) + Power(pz, 2))))) /
          (c_norm * p_dot_n_3 * tangential_dist);

      // d ["tangential cost"] / d [cx]
      jacobians[0][jacobian_index + 3] =
          (-2 * cx *
               (Power(cx - (c_dot_n * px) / p_dot_n, 2) +
                Power(cy - (c_dot_n * py) / p_dot_n, 2) +
                Power(cz - (c_dot_n * pz) / p_dot_n, 2)) +
           (2 * (Power(cx, 2) + Power(cy, 2) + Power(cz, 2)) *
            (-(cz * (-(nx * nz * Power(py, 2)) + Power(nx, 2) * px * pz +
                     Power(nz, 2) * px * pz + ny * py * (nz * px + nx * pz))) -
             cy * (Power(nx, 2) * px * py + nx * pz * (nz * py - ny * pz) +
                   ny * px * (ny * py + nz * pz)) +
             cx * (Power(ny * py + nz * pz, 2) +
                   Power(nx, 2) * (Power(py, 2) + Power(pz, 2))))) /
               Power(nx * px + ny * py + nz * pz, 2)) /
          (2.0 * Power(Power(cx, 2) + Power(cy, 2) + Power(cz, 2), 1.5) *
           tangential_dist);
      // d ["tangential cost"] / d [cy]
      jacobians[0][jacobian_index + 4] =
          (-2 * cy *
               (Power(cx - (c_dot_n * px) / p_dot_n, 2) +
                Power(cy - (c_dot_n * py) / p_dot_n, 2) +
                Power(cz - (c_dot_n * pz) / p_dot_n, 2)) -
           (2 * (Power(cx, 2) + Power(cy, 2) + Power(cz, 2)) *
            (cz * (Power(ny, 2) * py * pz + ny * px * (-(nz * px) + nx * pz) +
                   nz * py * (nx * px + nz * pz)) +
             cx * (Power(nx, 2) * px * py + nx * pz * (nz * py - ny * pz) +
                   ny * px * (ny * py + nz * pz)) -
             cy * (Power(nx, 2) * Power(px, 2) + 2 * nx * nz * px * pz +
                   Power(nz, 2) * Power(pz, 2) +
                   Power(ny, 2) * (Power(px, 2) + Power(pz, 2))))) /
               Power(nx * px + ny * py + nz * pz, 2)) /
          (2.0 * Power(Power(cx, 2) + Power(cy, 2) + Power(cz, 2), 1.5) *
           tangential_dist);
      // d ["tangential cost"] / d [cz]
      jacobians[0][jacobian_index + 5] =
          ((-2 * c_norm_2 *
            (-(cz * (Power(nx, 2) * Power(px, 2) + 2 * nx * ny * px * py +
                     Power(ny, 2) * Power(py, 2) +
                     Power(nz, 2) * (Power(px, 2) + Power(py, 2)))) +
             cx * (-(nx * nz * Power(py, 2)) + Power(nx, 2) * px * pz +
                   Power(nz, 2) * px * pz + ny * py * (nz * px + nx * pz)) +
             cy * (Power(ny, 2) * py * pz + ny * px * (-(nz * px) + nx * pz) +
                   nz * py * (nx * px + nz * pz)))) /
               Power(nx * px + ny * py + nz * pz, 2) -
           2 * cz *
               (Power(cx - (c_dot_n * px) / p_dot_n, 2) +
                Power(cy - (c_dot_n * py) / p_dot_n, 2) +
                Power(cz - (c_dot_n * pz) / p_dot_n, 2))) /
          (2.0 * Power(Power(cx, 2) + Power(cy, 2) + Power(cz, 2), 1.5) *
           tangential_dist);
      for (int i = 0; i < 6; ++i) {
        if (!std::isfinite(jacobians[0][jacobian_index + i])) {
          return false;
        }
      }

      if (!point_set_.weights.empty()) {
        float weight = has_weights ? std::sqrt(point_set_.weights[poi]) : 1.0f;
        for (int i = 0; i < 6; ++i) {
          // TODO(puneetl):  Optimize this!
          jacobians[0][jacobian_index + i] *= weight;
        }
      }
      jacobian_index += 6;
    }
  }
  return true;
}

}  // namespace tiler
}  // namespace seurat
