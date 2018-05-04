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

#ifndef VR_SEURAT_TILER_SELECTION_PARALLEL_SOLVER_H_
#define VR_SEURAT_TILER_SELECTION_PARALLEL_SOLVER_H_

#include <functional>
#include <memory>

#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_solver.h"

namespace seurat {
namespace tiler {
namespace selection {

// Given fixed multipliers, solves the resulting "dual" SelectionProblem by
// parallelizing over disjoint subexpressions.
//
// This is essentially a "dual decomposition" method.
class ParallelSolver : public DualSelectionSolver {
 public:
  typedef std::function<std::unique_ptr<DualSelectionSolver>(void)>
      SolverFactory;

  ParallelSolver(int thread_count, SolverFactory solver_factory);
  ~ParallelSolver() override;

  // DualSelectionSolver implementation.
  void Init(const SelectionProblem& problem) override;
  bool Solve(absl::Span<const double> multipliers,
             std::vector<int>* selected) override;

 private:
  // Scratch space for caching allocations.
  struct Workspace;
  const std::unique_ptr<Workspace> workspace_;

  // The maximum number of threads to use.
  const int thread_count_;

  // Constructs SelectionSolvers to use to solve problems derived from
  // subexpressions.
  const SolverFactory solver_factory_;

  // Solvers to solve subproblems with the given multipliers on different
  // threads.
  std::vector<std::unique_ptr<DualSelectionSolver>> dual_solvers_;

  // The current (global) problem being solved.
  SelectionProblem problem_;
};

}  // namespace selection
}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SELECTION_PARALLEL_SOLVER_H_
