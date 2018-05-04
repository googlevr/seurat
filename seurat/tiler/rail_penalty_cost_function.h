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

#ifndef VR_SEURAT_TILER_RAIL_PENALTY_COST_FUNCTION_H_
#define VR_SEURAT_TILER_RAIL_PENALTY_COST_FUNCTION_H_

#include <array>

#include "ceres/ceres.h"
#include "ion/math/range.h"
#include "ion/math/vector.h"

namespace seurat {
namespace tiler {

// Adds a quadratic penalty when the intersection of a plane and the given
// |rails| is outside the specified |ray_depth|.
//
// Put simply, this penalizes grazing-angle planes and enforces that resulting
// planes intersect the rails at finite values.
class RailPenaltyCostFunction : public ceres::CostFunction {
 public:
  RailPenaltyCostFunction(const std::array<ion::math::Vector3f, 4>& rails,
                          const ion::math::Range1f& ray_depth)
      : rails_(rails), ray_depth_(ray_depth) {}
  ~RailPenaltyCostFunction() override = default;

  // Perform initialization required to use this with ceres.
  //
  // If Evaluate() is being invoked outside of ceres, then invoking this is not
  // necessary.  Internally, this performs an allocation which can be very
  // expensive if performed for each call to Evaluate().
  void Initialize() {
    // 6 parameters.
    mutable_parameter_block_sizes()->push_back(6);
    // 2 residuals for each corner ray.
    set_num_residuals(8);
  }

  // CostFunction implementation.
  //
  // |parameters| must be:
  //   {{normal[0], normal[1], normal[2], center[0], center[1], center[2]}}
  // where the 'center' and 'normal' indicate a point and normal of a plane.
  //
  // A pair of residuals is returned for each corner ray with
  //   max(min_depth - t_hit, 0)
  //   min(t_hit - max_depth, 0)
  // where 't_hit' is the result of intersecting the corner ray with the
  // plane.
  bool Evaluate(double const* const* parameters, double* residuals,
                double** jacobians) const override;

 private:
  // Normalized rail directions.
  const std::array<ion::math::Vector3f, 4> rails_;

  // The range of ray depths in which no penalty is applied.
  const ion::math::Range1f ray_depth_;
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_RAIL_PENALTY_COST_FUNCTION_H_
