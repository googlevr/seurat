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

#include "seurat/tiler/selection/simplifying_solver.h"

#include <iostream>

#include "ion/base/logging.h"
#include "absl/types/span.h"
#include "seurat/base/util.h"
#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_util.h"

namespace seurat {
namespace tiler {
namespace selection {

void SimplifyingSolver::BuildSimplifiedProblem(
    absl::Span<const Token> original_expression,
    absl::Span<const int> subexpression_size, const ItemSet& original_items) {
  // Recursively translate |original_expression| and |original_items| into
  // |this->simplified_expression| and |this->simplified_items|.
  //
  // Specifically, this looks for clauses like And(a, b, c, ...) and simplifies
  // them into And(a', ...) where a' is a new "super-item" used to signify the
  // selection of a, b, and c.
  //
  // This saves time when the delegate solver must compute the partial-dual-cost
  // of selecting a, b, and c.
  //
  // Clauses which do not match this pattern are translated verbatim.

  if (original_expression.empty()) {
    // Base case.
    return;
  }
  // Process the subexpression starting from this token.
  //
  // Depending on whether this expression is a Token::Or(), Token::Item(), or
  // Token::And(), the translation is done differently.
  const Token& start = original_expression.front();
  CHECK_EQ(subexpression_size.front(), original_expression.size());
  if (start.IsOr()) {
    // Effectively do nothing here by simply outputting OR(...), recursing to
    // process subexpressions.
    simplified_expression_.push_back(Token::Or());
    int expression_offset = 1;
    while (expression_offset < original_expression.size() - 1) {
      int size = subexpression_size[expression_offset];
      BuildSimplifiedProblem(
          original_expression.subspan(expression_offset, size),
          subexpression_size.subspan(expression_offset, size), original_items);
      expression_offset += size;
    }
    simplified_expression_.push_back(Token::End());
  } else if (start.IsItem()) {
    // Effectively do nothing here by simply outputting Item(), but also link
    // together simplified_items_ to original_items.
    int old_item = start.GetItem();
    int new_item = simplified_items_.AppendItem(
        original_items.GetCost(old_item), original_items.GetWeights(old_item));
    simplified_expression_.push_back(Token::Item(new_item));
    old_items_from_simplified_item_.push_back({});
    CHECK_EQ(new_item, old_items_from_simplified_item_.size() - 1);
    old_items_from_simplified_item_.back().push_back(old_item);
  } else if (start.IsAnd()) {
    // For And() clauses, simplify as follows.
    //
    // If there are any immediate children of the And() clause which are items,
    // then batch them together as a single "super item" (with costs & weights
    // as the sum of all others).
    //
    // For all other children of the And() clause, recurse.
    simplified_expression_.push_back(Token::And());

    std::vector<int> child_items;

    // Iterate over subexpressions, recursing when appropriate.
    int expression_offset = 1;
    while (expression_offset < original_expression.size() - 1) {
      if (original_expression[expression_offset].IsItem()) {
        child_items.push_back(original_expression[expression_offset].GetItem());
        expression_offset++;
      } else {
        // If this isn't an item, simply recurse.
        int size = subexpression_size[expression_offset];
        BuildSimplifiedProblem(
            original_expression.subspan(expression_offset, size),
            subexpression_size.subspan(expression_offset, size),
            original_items);
        expression_offset += size;
      }
    }

    // Emit the "super item".
    if (!child_items.empty()) {
      float merged_cost = 0.0f;
      std::vector<float> merged_weights(original_items.GetNumWeights());
      for (int i : child_items) {
        merged_cost += original_items.GetCost(i);
        for (const auto& w : original_items.GetWeights(i)) {
          merged_weights[w.index] += w.value;
        }
      }
      std::vector<ItemSet::Weight> sparse_weights;
      for (int w = 0; w < merged_weights.size(); ++w) {
        if (merged_weights[w] > 0.0f) {
          ItemSet::Weight weight;
          weight.index = w;
          weight.value = merged_weights[w];
          sparse_weights.push_back(weight);
        }
      }
      int new_item = simplified_items_.AppendItem(merged_cost, sparse_weights);
      old_items_from_simplified_item_.push_back(std::move(child_items));
      CHECK_EQ(new_item, old_items_from_simplified_item_.size() - 1);
      simplified_expression_.push_back(Token::Item(new_item));
    }

    simplified_expression_.push_back(Token::End());
  }
}

void SimplifyingSolver::Init(const SelectionProblem& problem) {
  simplified_expression_.clear();
  simplified_items_ = ItemSet(problem.items->GetNumWeights());
  old_items_from_simplified_item_.clear();

  std::vector<int> sizes(problem.expression.size());
  PrecomputeSubexpressionSize(problem.expression, absl::MakeSpan(sizes));

  BuildSimplifiedProblem(problem.expression, sizes, *problem.items);

  SelectionProblem simplified;
  simplified.expression = simplified_expression_;
  simplified.capacity = problem.capacity;
  simplified.items = &simplified_items_;
  solver_->Init(simplified);
}

bool SimplifyingSolver::Solve(absl::Span<double> multipliers,
                              std::vector<int>* selected) {
  selected->clear();
  std::vector<int> simplified_solution;
  bool success = solver_->Solve(multipliers, &simplified_solution);
  if (!success) {
    return false;
  }
  for (int i : simplified_solution) {
    selected->insert(selected->end(),
                     old_items_from_simplified_item_[i].begin(),
                     old_items_from_simplified_item_[i].end());
  }
  return success;
}

}  // namespace selection
}  // namespace tiler
}  // namespace seurat
