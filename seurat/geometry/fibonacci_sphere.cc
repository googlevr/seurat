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

#include "seurat/geometry/fibonacci_sphere.h"

#include <math.h>  // For M_PI on MSVC14
#include <array>
#include <cmath>

#include "ion/math/matrix.h"
#include "ion/math/matrixutils.h"
#include "ion/math/vectorutils.h"

namespace seurat {
namespace geometry {

using ion::math::Matrix2d;
using ion::math::Point3d;
using ion::math::Point3f;
using ion::math::Vector2d;
using ion::math::Vector3d;

namespace {

double frac(double x) { return x - std::floor(x); }

const double kPhi = (1.0 + std::sqrt(5.0)) / 2.0;

// A fibonacci-sphere cell is the set of indices into a fibonacci-sphere point
// set representing the 4 closest points on the sphere to a given
// |normalized_direction| vector.
//
// For a nice visualization of this, see Figure 4 of the paper.
std::array<int, 4> FibonacciSphereCell(
    int num_points, double scrambler,
    const ion::math::Vector3d& normalized_direction) {
  constexpr double kPi = static_cast<double>(M_PI);
  double phi =
      std::min(std::atan2(normalized_direction[1], normalized_direction[0]),
               kPi) -
      scrambler;
  double cos_theta = normalized_direction[2];

  double k =
      std::max(2.0, std::floor(std::log(num_points * kPi * std::sqrt(5.0) *
                                        (1.0 - cos_theta * cos_theta)) /
                               std::log(kPhi * kPhi)));

  double fk = std::pow(kPhi, k) / std::sqrt(5.0);
  double f0 = std::round(fk);
  double f1 = std::round(fk * kPhi);

  Matrix2d b(2.0 * kPi * frac((f0 + 1.0) * (kPhi - 1.0)) -
                 2.0 * kPi * (kPhi - 1.0),  //
             2.0 * kPi * frac((f1 + 1.0) * (kPhi - 1.0)) -
                 2.0 * kPi * (kPhi - 1.0),  //
             -2.0 * f0 / num_points,        //
             -2.0 * f1 / num_points);

  Matrix2d inv_b = ion::math::Inverse(b);

  Vector2d c = inv_b * Vector2d(phi, cos_theta - (1.0 - 1.0 / num_points));
  c[0] = std::floor(c[0]);
  c[1] = std::floor(c[1]);

  std::array<int, 4> cell_indices;
  for (int s = 0; s < 4; ++s) {
    double cos_theta_unclamped =
        ion::math::Dot({b(1, 0), b(1, 1)}, Vector2d(s % 2, s / 2) + c) +
        (1.0 - 1.0 / num_points);
    double cos_theta = cos_theta_unclamped;
    if (cos_theta_unclamped > 1.0) {
      cos_theta = 2.0 - cos_theta;
    } else if (cos_theta_unclamped < -1.0) {
      cos_theta = -2.0 - cos_theta;
    } else {
      cos_theta = cos_theta * 2.0 - cos_theta;
    }

    int i = static_cast<int>(
        std::floor(num_points * 0.5 - cos_theta * num_points * 0.5));
    cell_indices[s] = i;
  }
  return cell_indices;
}

}  // namespace

Point3d GenerateFibonacciSpherePoint(int num_points, double scrambler, int i) {
  DCHECK_GT(num_points, 0);
  DCHECK_LT(i, num_points);
  double d_i = i;
  double phi = 2.0 * static_cast<double>(M_PI) * frac(d_i * (kPhi - 1.0));
  phi += scrambler;
  double cos_theta = 1.0 - (2.0 * d_i + 1) * (1.0 / num_points);
  double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
  return {std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, cos_theta};
}

int InverseFibonacciSphereMapping(
    int num_points, double scrambler,
    const ion::math::Vector3d& normalized_direction) {
  std::array<int, 4> cell_indices =
      FibonacciSphereCell(num_points, scrambler, normalized_direction);
  int best;
  double best_distance2 = std::numeric_limits<float>::infinity();
  for (const int cell_index : cell_indices) {
    double distance2 = ion::math::LengthSquared(
        GenerateFibonacciSpherePoint(num_points, scrambler, cell_index) -
        Point3d::Zero() - normalized_direction);
    if (distance2 < best_distance2) {
      best_distance2 = distance2;
      best = cell_index;
    }
  }
  return best;
}

}  // namespace geometry
}  // namespace seurat
