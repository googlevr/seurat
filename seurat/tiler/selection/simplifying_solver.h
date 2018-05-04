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

#ifndef VR_SEURAT_TILER_SELECTION_SIMPLIFYING_SOLVER_H_
#define VR_SEURAT_TILER_SELECTION_SIMPLIFYING_SOLVER_H_

#include <memory>
#include <vector>

#include "absl/types/span.h"
#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_solver.h"

namespace seurat {
namespace tiler {
namespace selection {

// Wraps another SelectionSolver, translating to and from a simplified
// SelectionProblem with fewer items and a smaller expression.
//
// Specifically, this looks for clauses like And(a, b, c, ...) and simplifies
// them into And(a', ...) where a' is a new "super-item" used to signify the
// selection of a, b, and c.
//
// This saves time when the delegate solver must compute the partial-dual-cost
// of selecting a, b, and c.
//
// Clauses which do not match this pattern are translated verbatim.
class SimplifyingSolver : public SelectionSolver {
 public:
  explicit SimplifyingSolver(std::unique_ptr<SelectionSolver> solver)
      : solver_(std::move(solver)) {}
  ~SimplifyingSolver() override = default;

  // SelectionSolver implementation.
  void Init(const SelectionProblem& problem) override;
  bool Solve(absl::Span<double> multipliers,
             std::vector<int>* selected) override;

 private:
  // Recursively builds the simplified problem, translating from the original
  // expression.
  void BuildSimplifiedProblem(absl::Span<const Token> original_expression,
                              absl::Span<const int> subexpression_size,
                              const ItemSet& original_items);

  // The solver to wrap.
  const std::unique_ptr<SelectionSolver> solver_;

  // The original problem being solved.
  SelectionProblem problem_;

  // A simplified version of the original problem.
  std::vector<Token> simplified_expression_;

  // The ItemSet corresponding to the simplified expression.
  ItemSet simplified_items_;

  // Returns the set of items in the original problem which correspond to each
  // item in the new, simplified, problem.
  std::vector<std::vector<int>> old_items_from_simplified_item_;
};

}  // namespace selection
}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SELECTION_SIMPLIFYING_SOLVER_H_
