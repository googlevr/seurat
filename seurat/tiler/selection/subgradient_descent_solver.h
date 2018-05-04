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

#ifndef VR_SEURAT_TILER_SELECTION_SUBGRADIENT_DESCENT_SOLVER_H_
#define VR_SEURAT_TILER_SELECTION_SUBGRADIENT_DESCENT_SOLVER_H_

#include <memory>

#include "absl/types/span.h"
#include "seurat/tiler/selection/cost_calculator.h"
#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_solver.h"

namespace seurat {
namespace tiler {
namespace selection {

// Solves SelectionProblems via subgradient-descent.
//
// Note that the final solution is *not* guaranteed to fit within the problem's
// capacity. In other words, it may not be primal-feasible.
//
// The implementation follows directly from "An Applications Oriented Guide
// to Lagrangian Relaxation" (Fisher, 1985) and "The Lagrangian Relaxation
// Method for Solving Integer Programming Problems" (Fisher, 2004).
class SubgradientDescentSolver : public SelectionSolver {
 public:
  SubgradientDescentSolver(int thread_count,
                           std::unique_ptr<DualSelectionSolver> dual_solver,
                           int iteration_count);
  ~SubgradientDescentSolver() override;

  // SelectionSolver implementation.
  void Init(const SelectionProblem& problem) override;
  bool Solve(absl::Span<double> multipliers,
             std::vector<int>* selected) override;

 private:
  // Scratch space for caching allocations.
  struct Workspace;
  const std::unique_ptr<Workspace> workspace_;

  // A solver to solve subproblems with fixed multipliers.
  const std::unique_ptr<DualSelectionSolver> dual_solver_;

  // The number of iterations of subgradient-descent to perform.
  const int iteration_count_;

  // The number of iterations without improvement before the subgradient-descent
  // step-size is halved.
  const int num_iterations_before_halving_step_size_;

  // The current problem being solved.
  SelectionProblem problem_;

  // Used to compute total costs & weight of intermediate solutions.
  CostCalculator cost_calculator_;
};

}  // namespace selection
}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SELECTION_SUBGRADIENT_DESCENT_SOLVER_H_
