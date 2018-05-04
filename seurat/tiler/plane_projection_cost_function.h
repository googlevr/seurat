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

#ifndef VR_SEURAT_TILER_PLANE_PROJECTION_COST_FUNCTION_H_
#define VR_SEURAT_TILER_PLANE_PROJECTION_COST_FUNCTION_H_

#include "ceres/ceres.h"
#include "absl/types/span.h"
#include "seurat/tiler/point_set.h"

namespace seurat {
namespace tiler {

// Models the distance between projections of points onto planes.
//
// See Evaluate() for details.
class PlaneProjectionCostFunction : public ceres::CostFunction {
 public:
  PlaneProjectionCostFunction(absl::Span<const int> points_of_interest,
                              const PointSet& point_set)
      : points_of_interest_(points_of_interest), point_set_(point_set) {}
  ~PlaneProjectionCostFunction() override = default;

  // Perform initialization required to use this with ceres.
  //
  // If Evaluate() is being invoked outside of ceres, then invoking this is not
  // necessary.  Internally, this performs an allocation which can be very
  // expensive if performed for each call to Evaluate().
  void Initialize() {
    mutable_parameter_block_sizes()->push_back(6);
    // 1 residual for each point of interest.
    set_num_residuals(points_of_interest_.size());
  }

  // CostFunction implementation.
  //
  // |parameters| must be:
  //   {{normal[0], normal[1], normal[2], center[0], center[1], center[2]}}
  // where the 'center' and 'normal' indicate a point and normal of a plane.
  //
  // A residual is returned for each point-of-interest with
  //   (1 - t_hit) * sqrt(weight)
  // where 't_hit' is the result of intersecting the Origin->Point ray with the
  // plane and 'weight' is the weight of the point.
  bool Evaluate(double const* const* parameters, double* residuals,
                double** jacobians) const override;

 private:
  // The points in the PointSet to which we are fitting the proxy geometry.
  const absl::Span<const int> points_of_interest_;

  // The set of all points.
  const PointSet point_set_;
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_PLANE_PROJECTION_COST_FUNCTION_H_
