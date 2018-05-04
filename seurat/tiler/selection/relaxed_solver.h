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

#ifndef VR_SEURAT_TILER_SELECTION_RELAXED_SOLVER_H_
#define VR_SEURAT_TILER_SELECTION_RELAXED_SOLVER_H_

#include <memory>

#include "absl/types/span.h"
#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_solver.h"

namespace seurat {
namespace tiler {
namespace selection {

// An optimized, single-threaded DualSelectionSolver.
//
// Like all DualSelectionSolvers, this optimizes the "relaxed" problem,
// evaluated at the given multipliers. See selection_problem.h for details.
class RelaxedSolver : public DualSelectionSolver {
 public:
  RelaxedSolver();
  ~RelaxedSolver() override;

  // DualSelectionSolver implementation.
  void Init(const SelectionProblem& problem) override;
  bool Solve(absl::Span<const double> multipliers,
             std::vector<int>* selected) override;

 private:
  // Temporary allocations needed for Solve().
  struct Workspace;
  std::unique_ptr<Workspace> workspace_;

  // The current problem to be solved.
  SelectionProblem problem_;
};

}  // namespace selection
}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SELECTION_RELAXED_SOLVER_H_
