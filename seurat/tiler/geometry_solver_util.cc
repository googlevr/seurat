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

#include "seurat/tiler/geometry_solver_util.h"

#include "absl/types/span.h"
#include "seurat/base/parallel.h"

namespace seurat {
namespace tiler {

using ion::math::Point3f;
using ion::math::Vector3f;

Point3f ComputeInitialCenterPoint(const PointSet& point_set,
                                  absl::Span<const int> points_of_interest) {
  Vector3f weighted_mean_pos;
  float total_weight = 0.0f;
  for (const int i : points_of_interest) {
    const Point3f& point = point_set.positions[i];
    float weight = point_set.weights.empty() ? 1.0f : point_set.weights[i];
    weight /= ion::math::Length(point - Point3f::Zero());
    weighted_mean_pos = weighted_mean_pos + (point * weight - Point3f::Zero());
    total_weight += weight;
  }
  weighted_mean_pos /= total_weight;
  return weighted_mean_pos + Point3f::Zero();
}

bool InitialEvaluationSucceeds(ceres::Problem* problem) {
  double cost;
  std::vector<double> residuals;
  std::vector<double> gradient;
  // Request residuals & gradients to ensure the full Evaluate() function is
  // invoked.
  return problem->Evaluate(ceres::Problem::EvaluateOptions(), &cost, &residuals,
                           &gradient, /*jacobian=*/nullptr);
}

}  // namespace tiler
}  // namespace seurat
