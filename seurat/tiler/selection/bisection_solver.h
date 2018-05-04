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

#ifndef VR_SEURAT_TILER_SELECTION_BISECTION_SOLVER_H_
#define VR_SEURAT_TILER_SELECTION_BISECTION_SOLVER_H_

#include <memory>

#include "absl/types/span.h"
#include "seurat/tiler/selection/cost_calculator.h"
#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_solver.h"

namespace seurat {
namespace tiler {
namespace selection {

// Use bisection to evenly scale all multipliers until the primal is feasible.
//
// Note that this is *not* a standard technique.  However, it terminates
// in a fixed number of iterations and works well in practice. Fisher (2004)
// calls this kind of technique a "Lagrangian heuristic" method.
class BisectionSolver : public SelectionSolver {
 public:
  BisectionSolver(int thread_count, std::unique_ptr<DualSelectionSolver> solver,
                  int iteration_count, bool initialize_multipliers = false);
  ~BisectionSolver() override;

  // SelectionSolver implementation.
  void Init(const SelectionProblem& problem) override;
  bool Solve(absl::Span<double> multipliers,
             std::vector<int>* selected) override;

 private:
  // Scratch space for caching allocations.
  struct Workspace;
  const std::unique_ptr<Workspace> workspace_;

  // A solver to solve subproblems with fixed multipliers.
  const std::unique_ptr<DualSelectionSolver> solver_;

  // The number of iterations of bisection to perform.
  const int iteration_count_;

  // Computes cost & weight for solutions.
  CostCalculator cost_calculator_;

  // The current problem being solved.
  SelectionProblem problem_;

  // Whether this solver should ignore existing multipliers and initialize them.
  bool initialize_multipliers_;
};

}  // namespace selection
}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SELECTION_BISECTION_SOLVER_H_
