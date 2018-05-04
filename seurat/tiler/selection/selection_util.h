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

#ifndef VR_SEURAT_TILER_SELECTION_SELECTION_UTIL_H_
#define VR_SEURAT_TILER_SELECTION_SELECTION_UTIL_H_

#include "seurat/tiler/selection/selection_problem.h"

namespace seurat {
namespace tiler {
namespace selection {

// Returns true if |result| has a weight within the |problem|'s capacity.
bool IsFeasibleWeight(const SelectionProblem& problem,
                      const absl::Span<const double> total_weight);

// Computes the total cost of all |items| from the given |item_set|.
double ComputeTotalCost(const ItemSet& item_set, absl::Span<const int> items);

// Computes the total weight of all |items| and returns it as a dense vector in
// |weight|.
void ComputeTotalWeight(const ItemSet& item_set, absl::Span<const int> items,
                        absl::Span<double> weight);

// Computes the total cost of the dual to the SelectionProblem, evaluated with
// the given multipliers, assuming the |weight| and |cost| are the result of
// selecting the optimal items for these |weight_multipliers|.
double ComputeDualCost(const SelectionProblem& problem,
                       absl::Span<const double> weight_multipliers,
                       absl::Span<const double> weight, double cost);

// Returns the given |expression| as a human-readable string.
std::string PrintExpressionToString(absl::Span<const Token> expression);

// Parses the |expression| to compute the size (the number of tokens) of each
// subexpression.
//
// For example, the expression
//   AND(1, 2, 3, OR(4))
// corresponds to the sequence:
//   AND, 1, 2, 3, OR, 4, END, END
// and has subexpression sizes:
//   8, 1, 1, 1, 3, 1, 3, 8
//
// For example, the first token, AND, is a clause containing 8 total tokens.
void PrecomputeSubexpressionSize(absl::Span<const Token> expression,
                                 absl::Span<int> subexpression_size);

// Returns whether |selected_items| are a valid selection, satisfying the
// constraints implied by the |expression|.
//
// Note that empty expressions, such as AND() and OR(), are never valid.
bool ValidateSelection(absl::Span<const Token> expression,
                       absl::Span<const int> selected_items);

}  // namespace selection
}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SELECTION_SELECTION_UTIL_H_
