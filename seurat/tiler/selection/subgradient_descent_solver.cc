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

#include "seurat/tiler/selection/subgradient_descent_solver.h"

#include <cmath>
#include <cstdlib>

#include "absl/types/span.h"
#include "seurat/base/util.h"
#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_util.h"

namespace seurat {
namespace tiler {
namespace selection {

struct SubgradientDescentSolver::Workspace {
  std::vector<double> best_multipliers;
  std::vector<double> current_multipliers;
  std::vector<double> current_weight;
};

SubgradientDescentSolver::SubgradientDescentSolver(
    int thread_count, std::unique_ptr<DualSelectionSolver> dual_solver,
    int iteration_count)
    : workspace_(new Workspace),
      dual_solver_(std::move(dual_solver)),
      iteration_count_(iteration_count),
      num_iterations_before_halving_step_size_(5),
      cost_calculator_(thread_count) {}
SubgradientDescentSolver::~SubgradientDescentSolver() {}

void SubgradientDescentSolver::Init(const SelectionProblem& problem) {
  problem_ = problem;
  dual_solver_->Init(problem);
}

bool SubgradientDescentSolver::Solve(absl::Span<double> multipliers,
                                     std::vector<int>* selected) {
  const int num_weights = problem_.items->GetNumWeights();
  selected->clear();

  workspace_->best_multipliers.resize(num_weights);
  workspace_->current_multipliers.resize(num_weights);
  workspace_->current_weight.resize(num_weights);

  // Pointers to the multipliers being evaluated and the best multipliers found
  // so far.
  //
  // As better values are found, these pointers will swap back and forth.
  std::vector<double>* current_multipliers = &workspace_->current_multipliers;
  std::vector<double>* best_multipliers = &workspace_->best_multipliers;

  for (int d = 0; d < num_weights; ++d) {
    current_multipliers->at(d) = multipliers[d];
  }

  // This implementation follows *directly* from the literature, based on the
  // classic subgradient-descent method.
  //
  // See Fisher (1985) or Fisher (2004) for details.

  // The cost of the dual for the best value found so far.
  //
  // This is what we're trying to improve with subgradient descent.
  //
  // The multipliers at which this was computed are in best_multipliers.
  double best_dual_cost = -std::numeric_limits<double>::infinity();

  // Initialize the best multipliers with the initial values.
  std::copy(current_multipliers->begin(), current_multipliers->end(),
            best_multipliers->begin());

  // The iteration at which we achieved the best_solution so far, used to
  // adaptively scale the step size.
  int best_solution_iteration = -1;

  // The (primal) cost of the best feasible solution found so far.
  double best_feasible_primal_cost = -1.0;

  // The step-size scale-factor (between 0.0 and 2.0).
  double alpha = 2.0;

  for (int iter = 0; iter < iteration_count_; ++iter) {
    // The multipliers being evaluated in this iteration.
    absl::Span<const double> evaluated_multipliers = *current_multipliers;
    selected->clear();
    if (!dual_solver_->Solve(evaluated_multipliers, selected)) {
      return false;
    }

    double current_dual_cost;
    double current_primal_cost;
    cost_calculator_.ComputeSolutionCost(
        problem_, evaluated_multipliers, *selected, &current_primal_cost,
        &current_dual_cost, absl::MakeSpan(workspace_->current_weight));

    if (IsFeasibleWeight(problem_, workspace_->current_weight)) {
      best_feasible_primal_cost =
          std::min(best_feasible_primal_cost, current_primal_cost);
    }

    if (current_dual_cost > best_dual_cost) {
      best_dual_cost = current_dual_cost;

      std::vector<double>* tmp = best_multipliers;
      best_multipliers = current_multipliers;
      current_multipliers = tmp;

      best_solution_iteration = iter;
    }

    // Adjust the step-size according to the classic adaptive method attributed
    // to Held et al. (1974).
    if (iter - best_solution_iteration >
        num_iterations_before_halving_step_size_) {
      alpha /= 2.0;
      // Reset the iteration-counter to ensure we do not keep halving the
      // step size if the next step fails to improve the dual.
      best_solution_iteration = iter;
    }
    double step_size_denom = 0.0;
    for (int d = 0; d < num_weights; ++d) {
      double slack = workspace_->current_weight[d] - problem_.capacity[d];
      step_size_denom += slack * slack;
    }
    double step_size =
        alpha * std::fabs(best_feasible_primal_cost - current_dual_cost) /
        step_size_denom;

    // Adjust the multipliers using the standard subgradient-descent method.
    for (int d = 0; d < num_weights; ++d) {
      current_multipliers->at(d) =
          std::max(0.0, evaluated_multipliers[d] +
                            step_size * (workspace_->current_weight[d] -
                                         problem_.capacity[d]));
    }
  }

  std::copy(best_multipliers->begin(), best_multipliers->end(),
            multipliers.begin());
  return dual_solver_->Solve(multipliers, selected);
}

}  // namespace selection
}  // namespace tiler
}  // namespace seurat
