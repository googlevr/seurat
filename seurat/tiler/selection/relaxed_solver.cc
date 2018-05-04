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

#include "seurat/tiler/selection/relaxed_solver.h"

#include <iostream>
#include <numeric>
#include <vector>

#include "ion/base/logging.h"
#include "absl/types/span.h"
#include "seurat/base/util.h"
#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_util.h"

namespace seurat {
namespace tiler {
namespace selection {

namespace {

// The partial-dual-cost is the component of the dual which is a function of the
// given variable.
double ComputePartialDualCost(const ItemSet& items,
                              absl::Span<const double> multipliers, int item) {
  double cost = items.GetCost(item);
  for (const auto& w : items.GetWeights(item)) {
    cost += multipliers[w.index] * w.value;
  }
  return cost;
}

struct PartialSolution {
  int item_count;
  double partial_dual_cost;
  bool success;
};

}  // namespace

struct RelaxedSolver::Workspace {
  std::vector<int> subexpression_size;
  std::vector<Token> path;
  std::vector<int> path_end;
  std::vector<double> max_partial_dual_cost;
  std::vector<PartialSolution> subtree_solution;
};

RelaxedSolver::RelaxedSolver() : workspace_(new Workspace) {}
RelaxedSolver::~RelaxedSolver() {}

void RelaxedSolver::Init(const SelectionProblem& problem) {
  problem_ = problem;
  workspace_->subexpression_size.resize(problem.expression.size());
  PrecomputeSubexpressionSize(problem.expression,
                              absl::MakeSpan(workspace_->subexpression_size));
}

bool RelaxedSolver::Solve(absl::Span<const double> multipliers,
                          std::vector<int>* selected) {
  // This is the inner-loop of the SelectionSolver.
  //
  // Optimizing for performance here is critical.
  selected->clear();

  absl::Span<const Token> expression = problem_.expression;
  absl::Span<const int> subexpression_size = workspace_->subexpression_size;
  const ItemSet& items = *problem_.items;

  if (expression.empty()) return false;

  // This method proceeds as follows:
  // Process each token one at a time.
  // Maintain, on a logical stack, (partial) solutions to subtrees:
  //   * If the parent node is AND, then the subtree solution consists of the
  //     merged solutions to all of its child nodes.
  //   * If the parent node is OR, then the subtree solution consists of the
  //     current best solution from all child nodes processed so far.

  // A stack maintaining the current parent node, either And() or Or().
  std::vector<Token>& path = workspace_->path;
  path.clear();
  // The index of the End() corresponding to the node at path.back().
  std::vector<int>& path_end = workspace_->path_end;
  path_end.clear();
  // The maximum partial-dual-cost of any solution to the subtree defined by
  // path.back().
  std::vector<double>& max_partial_dual_cost =
      workspace_->max_partial_dual_cost;
  max_partial_dual_cost.clear();
  max_partial_dual_cost.push_back(std::numeric_limits<double>::infinity());

  // For each subtree, this is the current number of items which have been
  // appended to |selected| for that subtree's solution.
  std::vector<PartialSolution>& subtree_solution = workspace_->subtree_solution;
  subtree_solution.clear();

  // Loop until the last End() token.
  for (int expression_index = 0; expression_index < expression.size() - 1;
       ++expression_index) {
    const Token& token = expression[expression_index];
    if (token.IsAnd() || token.IsOr()) {
      path.push_back(token);
      path_end.push_back(expression_index +
                         subexpression_size[expression_index] - 1);
      max_partial_dual_cost.push_back(max_partial_dual_cost.back());
      if (token.IsAnd()) {
        // Start with a partial-dual-cost of 0, since And(a, b, c) is the sum of
        // the partial-dual-costs of a, b, and c.
        //
        // Start with success, since the success of And(a, b, c) is the
        // logical-and of the success of a, b, and c. (Unless the And() clause
        // is empty, in which case it is marked as failure below!)
        subtree_solution.push_back({0, 0.0, true});
      } else {
        // Start with a partial-dual-cost of Infinity, since Or(a, b, c) is the
        // min of the partial-dual-costs of a, b, and c.
        //
        // Start with failure, since the success of Or(a, b, c) is the
        // logical-or of the success of a, b, and c.
        subtree_solution.push_back(
            {0, std::numeric_limits<double>::infinity(), false});
      }
    } else if (token.IsItem()) {
      if (path.back().IsAnd()) {
        // Add this item to the current solution.
        const int item = token.GetItem();
        selected->push_back(item);
        subtree_solution.back().item_count++;
        subtree_solution.back().partial_dual_cost +=
            ComputePartialDualCost(items, multipliers, item);
        if (subtree_solution.back().partial_dual_cost >
            max_partial_dual_cost.back()) {
          // The current solution is now so expensive that it'll never be
          // selected, so early-exit to the end of this And() clause by
          // skipping ahead in the expression.
          //
          // This early-exit saves ~30% runtime on real problems.
          expression_index = path_end.back() - 1;
          subtree_solution.back().success = false;
        }
      } else {
        CHECK(path.back().IsOr());
        const int item = token.GetItem();
        double item_cost = ComputePartialDualCost(items, multipliers, item);
        if (item_cost <= subtree_solution.back().partial_dual_cost) {
          // Replace the existing solution with this one.
          selected->resize(selected->size() -
                           subtree_solution.back().item_count);
          selected->push_back(item);
          subtree_solution.back().item_count = 1;
          subtree_solution.back().partial_dual_cost = item_cost;
          subtree_solution.back().success = true;
          max_partial_dual_cost.back() =
              std::min(max_partial_dual_cost.back(), item_cost);
        }
      }
    } else if (token.IsEnd()) {
      if (path.back().IsAnd() && subtree_solution.back().item_count == 0) {
        // This is an empty And() clause, so mark failure.
        subtree_solution.back().success = false;
      }
      // Merge subtree_solution.back() into the one just before it, based on its
      // type.
      if (path[path.size() - 2].IsAnd()) {
        PartialSolution& parent_solution =
            subtree_solution[subtree_solution.size() - 2];
        parent_solution.item_count += subtree_solution.back().item_count;
        parent_solution.partial_dual_cost +=
            subtree_solution.back().partial_dual_cost;
        parent_solution.success &= subtree_solution.back().success;

        path.pop_back();
        path_end.pop_back();
        subtree_solution.pop_back();
        max_partial_dual_cost.pop_back();

        if (subtree_solution.back().partial_dual_cost >
            max_partial_dual_cost.back()) {
          // The current solution is now so expensive that it'll never be
          // selected, so early-exit to the end of this And() clause by
          // skipping ahead in the expression.
          expression_index = path_end.back() - 1;
          subtree_solution.back().success = false;
        }
      } else {
        CHECK(path[path.size() - 2].IsOr());
        PartialSolution& parent_solution =
            subtree_solution[subtree_solution.size() - 2];
        if (subtree_solution.back().success &&
            (!parent_solution.success ||
             subtree_solution.back().partial_dual_cost <
                 parent_solution.partial_dual_cost)) {
          // Replace the items from the existing solution with this one by
          // shifting this solution's items left to overwrite the old items.
          const int items_in_old_solution = parent_solution.item_count;
          const int items_in_this_solution = subtree_solution.back().item_count;
          const int items_begin = selected->size() - (items_in_this_solution +
                                                      items_in_old_solution);
          std::copy(selected->begin() + items_begin + items_in_old_solution,
                    selected->end(), selected->begin() + items_begin);
          selected->resize(selected->size() - items_in_old_solution);

          parent_solution.item_count = items_in_this_solution;
          parent_solution.partial_dual_cost =
              subtree_solution.back().partial_dual_cost;
          parent_solution.success = subtree_solution.back().success;
        } else {
          // Remove the current solution.
          selected->resize(selected->size() -
                           subtree_solution.back().item_count);
        }
        path.pop_back();
        path_end.pop_back();
        subtree_solution.pop_back();
        max_partial_dual_cost.pop_back();
      }
    }
  }
  CHECK_EQ(max_partial_dual_cost.size(), 2);
  CHECK_EQ(path.size(), 1);
  CHECK_EQ(path_end.size(), 1);
  CHECK_EQ(subtree_solution.size(), 1);
  CHECK_EQ(selected->size(), subtree_solution.back().item_count);

  bool success = subtree_solution.back().success && !selected->empty();

  DCHECK(!success || ValidateSelection(expression, *selected));

  return success;
}

}  // namespace selection
}  // namespace tiler
}  // namespace seurat
