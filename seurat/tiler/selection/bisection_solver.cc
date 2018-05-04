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

#include "seurat/tiler/selection/bisection_solver.h"

#include <cmath>
#include <numeric>

#include "absl/types/span.h"
#include "seurat/base/util.h"
#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_util.h"

namespace seurat {
namespace tiler {
namespace selection {

// Various dynamically-allocated arrays used by Solve().
//
// Stored here to cache allocations.
struct BisectionSolver::Workspace {
  std::vector<double> current_multipliers;
  std::vector<double> total_weight;
};

BisectionSolver::BisectionSolver(int thread_count,
                                 std::unique_ptr<DualSelectionSolver> solver,
                                 int iteration_count,
                                 bool initialize_multipliers)
    : workspace_(new Workspace),
      solver_(std::move(solver)),
      iteration_count_(iteration_count),
      cost_calculator_(thread_count),
      initialize_multipliers_(initialize_multipliers) {}
BisectionSolver::~BisectionSolver() {}

void BisectionSolver::Init(const SelectionProblem& problem) {
  problem_ = problem;
  solver_->Init(problem);
}

bool BisectionSolver::Solve(absl::Span<double> multipliers,
                            std::vector<int>* selected) {
  const int num_weights = problem_.capacity.size();
  selected->clear();

  // Scratch space to store the dense vector of multipliers currently being
  // evaluated.
  workspace_->current_multipliers.resize(num_weights);
  absl::Span<double> current_multipliers =
      absl::MakeSpan(workspace_->current_multipliers);

  // Scratch space to store the dense vector of the total weight of the most
  // recent solution.
  workspace_->total_weight.resize(num_weights);
  absl::Span<double> total_weight = absl::MakeSpan(workspace_->total_weight);

  if (initialize_multipliers_) {
    for (int d = 0; d < num_weights; ++d) {
      multipliers[d] = 1.0 / problem_.capacity[d];
    }
  }

  double lb = 0.0;
  double ub = 1.0;

  // Double the upper-bound until it results in a feasible solution.
  //
  // Probing in this way ensures we start the bisection step (below) with a
  // numerically-stable value, if possible.
  for (int i = 0; i < iteration_count_; ++i) {
    for (int d = 0; d < problem_.items->GetNumWeights(); ++d) {
      current_multipliers[d] = ub * multipliers[d];
      if (!std::isfinite(current_multipliers[d])) {
        // If we get here, everything is scaled so much that there's no hope of
        // finding a feasible set of multipliers.
        return false;
      }
    }
    if (!solver_->Solve(current_multipliers, selected)) {
      return false;
    }
    double primal_cost;
    double dual_cost;
    cost_calculator_.ComputeSolutionCost(problem_, current_multipliers,
                                         *selected, &primal_cost, &dual_cost,
                                         total_weight);
    if (IsFeasibleWeight(problem_, total_weight)) {
      break;
    }
    ub *= 2.0;
  }

  if (!IsFeasibleWeight(problem_, total_weight)) {
    return false;
  }

  // Bisect to find the best feasible upper bound.
  for (int i = 0; i < iteration_count_; ++i) {
    double pivot = (lb + ub) / 2.0;

    for (int d = 0; d < problem_.items->GetNumWeights(); ++d) {
      current_multipliers[d] = pivot * multipliers[d];
    }

    if (!solver_->Solve(current_multipliers, selected)) {
      return false;
    }
    double primal_cost;
    double dual_cost;
    cost_calculator_.ComputeSolutionCost(problem_, current_multipliers,
                                         *selected, &primal_cost, &dual_cost,
                                         total_weight);
    if (IsFeasibleWeight(problem_, total_weight)) {
      ub = pivot;
    } else {
      lb = pivot;
    }
  }

  // Run the solver using the feasible upper-bound we found.
  for (int d = 0; d < problem_.items->GetNumWeights(); ++d) {
    multipliers[d] = ub * multipliers[d];
  }
  return solver_->Solve(multipliers, selected);
}

}  // namespace selection
}  // namespace tiler
}  // namespace seurat
