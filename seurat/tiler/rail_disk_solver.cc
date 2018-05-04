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

#include "seurat/tiler/rail_disk_solver.h"

#include <array>
#include <cmath>
#include <vector>

#include "ceres/ceres.h"
#include "ion/math/vector.h"
#include "seurat/tiler/geometry_solver_util.h"
#include "seurat/tiler/plane_projection_cost_function.h"
#include "seurat/tiler/rail_penalty_cost_function.h"
#include "seurat/tiler/tangential_disk_cost_function.h"

namespace seurat {
namespace tiler {

using ion::math::Point3f;
using ion::math::Vector3f;

void RailDiskSolver::Init(const PointSet& point_set) {
  point_set_ = point_set;
  subdivision_->Init(point_set);
}

float RailDiskSolver::ComputeError(int point_index,
                                   const GeometryModel& model) const {
  // Add up the plane-projection & tangential costs.
  std::array<double, 6> params = {{
      model.normal[0], model.normal[1], model.normal[2],  //
      model.center[0], model.center[1], model.center[2]   //
  }};
  std::array<double const*, 1> parameter_block = {{params.data()}};
  std::array<double, 1> residuals;

  double total_cost = 0.0;
  PlaneProjectionCostFunction(absl::Span<const int>(&point_index, 1),
                              point_set_)
      .Evaluate(parameter_block.data(), residuals.data(), nullptr);
  total_cost += residuals[0] * residuals[0];

  if (tangential_factor_ != 0.0f) {
    TangentialDiskCostFunction(absl::Span<const int>(&point_index, 1),
                               point_set_)
        .Evaluate(parameter_block.data(), residuals.data(), nullptr);
    // TODO(puneetl):  That the tangential_factor_ is squared and multiplied by
    // 2 is largely to maintain compatibility with the original DiskSolver.
    // However, since it is a constant, it could be folded into a single value.
    total_cost += residuals[0] * residuals[0] * 2.0 *
                  (tangential_factor_ * tangential_factor_);
  }

  // The RailPenaltyCostFunction is intentionally absent, since it is
  // not a function of indvidual points.
  //
  // It is only applied in FitModel().

  return total_cost;
}

void RailDiskSolver::InitializeModel(int point_index,
                                     GeometryModel* model) const {
  model->center = point_set_.positions[point_index];
  model->normal = ion::math::Normalized(model->center - Point3f::Zero());
}

bool RailDiskSolver::FitModel(absl::Span<const int> point_indices,
                              GeometryModel* model) const {
  if (point_indices.size() < 3) {
    return false;
  }

  Point3f weighted_mean_pos =
      ComputeInitialCenterPoint(point_set_, point_indices);

  // The parameters to solve for.  As a starting point, use the
  // previously-estimated normal and a weighted average of the point positions.
  std::array<double, 6> params = {{
      model->normal[0],      //
      model->normal[1],      //
      model->normal[2],      //
      weighted_mean_pos[0],  //
      weighted_mean_pos[1],  //
      weighted_mean_pos[2]   //
  }};

  ceres::Problem problem;

  std::unique_ptr<PlaneProjectionCostFunction> projection_cost(
      new PlaneProjectionCostFunction(point_indices, point_set_));
  projection_cost->Initialize();
  problem.AddResidualBlock(projection_cost.release(), nullptr, params.data());

  if (tangential_factor_ != 0.0f) {
    std::unique_ptr<TangentialDiskCostFunction> tangential_cost(
        new TangentialDiskCostFunction(point_indices, point_set_));
    tangential_cost->Initialize();
    problem.AddResidualBlock(
        tangential_cost.release(),
        new ceres::ScaledLoss(nullptr,
                              tangential_factor_ * tangential_factor_ * 2.0,
                              ceres::DO_NOT_TAKE_OWNERSHIP),
        params.data());
  }

  std::array<Vector3f, 4> rails = subdivision_->GetRails(model->cell);
  double scale = point_indices.size();
  std::unique_ptr<RailPenaltyCostFunction> corner_ray_penalty(
      new RailPenaltyCostFunction(rails, depth_range_));
  corner_ray_penalty->Initialize();
  problem.AddResidualBlock(
      corner_ray_penalty.release(),
      new ceres::ScaledLoss(nullptr, scale, ceres::DO_NOT_TAKE_OWNERSHIP),
      params.data());

  if (!InitialEvaluationSucceeds(&problem)) return false;

  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_QR;
  options.num_threads = 1;
  options.logging_type = ceres::LoggingType::SILENT;

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  if (!summary.IsSolutionUsable()) {
    return false;
  }

  const Vector3f final_normal(params[0], params[1], params[2]);
  const Point3f final_center(params[3], params[4], params[5]);
  model->center = final_center;
  model->normal = ion::math::Normalized(final_normal);
  return true;
}

}  // namespace tiler
}  // namespace seurat
