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

#ifndef VR_SEURAT_TILER_SELECTION_SEQUENTIAL_SOLVER_H_
#define VR_SEURAT_TILER_SELECTION_SEQUENTIAL_SOLVER_H_

#include <memory>
#include <vector>

#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_solver.h"

namespace seurat {
namespace tiler {
namespace selection {

// Invokes a sequence of SelectionSolvers which may iteratively improve upon the
// multipliers, returning the solution from the last stage.
class SequentialSolver : public SelectionSolver {
 public:
  explicit SequentialSolver(
      std::vector<std::unique_ptr<SelectionSolver>> sequence)
      : sequence_(std::move(sequence)) {}
  ~SequentialSolver() override = default;

  // SelectionSolver implementation.
  void Init(const SelectionProblem& problem) override {
    for (auto& solver : sequence_) {
      solver->Init(problem);
    }
  }
  bool Solve(absl::Span<double> multipliers,
             std::vector<int>* selected) override {
    // Initial solvers in the sequence may fail.  This is okay so long as the
    // last result is a success.
    bool last_result = false;
    for (auto& solver : sequence_) {
      selected->clear();
      last_result = solver->Solve(multipliers, selected);
    }
    return last_result;
  }

 private:
  // The sequence of solvers to run.
  std::vector<std::unique_ptr<SelectionSolver>> sequence_;
};

}  // namespace selection
}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SELECTION_SEQUENTIAL_SOLVER_H_
