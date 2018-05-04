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

#include "seurat/tiler/selection/parallel_solver.h"

#include <iostream>
#include <queue>

#include "ion/base/logging.h"
#include "absl/types/span.h"
#include "seurat/base/parallel.h"
#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_util.h"

namespace seurat {
namespace tiler {
namespace selection {

namespace {

// Repeatedly splits the |expression| into non-overlapping subexpressions.
//
// Returns these subexpressions, which are slices into the input |expression|
// consisting of And() or Or() clauses, sorted according to their position in
// the expression.
//
// Requires a precomputed table of |subexpression_size|.  See
// PrecomputeSubexpressionSize() for details.
std::vector<absl::Span<const Token>> SplitIntoSubexpressions(
    absl::Span<const Token> expression,
    absl::Span<const int> subexpression_size, int min_expression_count) {
  struct Subexpression {
    absl::Span<const Token> slice;

    // Order by size.
    bool operator<(const Subexpression& rhs) const {
      return slice.size() < rhs.slice.size();
    }
  };

  // Repeatedly split the largest subexpression into nonoverlapping pieces.
  std::priority_queue<Subexpression> subexpressions;
  subexpressions.push(Subexpression{expression});
  while (subexpressions.size() < min_expression_count) {
    absl::Span<const Token> to_split = subexpressions.top().slice;
    subexpressions.pop();
    bool has_children = false;
    // Loop over all immediate children.
    //
    // Ignore the first and last token which must correspond to And()/Or() and
    // End(), respectively.
    int offset = 1;
    while (offset < to_split.size() - 1) {
      // The size of the expression starting at 'offset'.
      int size = subexpression_size[&(to_split[offset]) - expression.begin()];
      if (to_split[offset].IsAnd() || to_split[offset].IsOr()) {
        has_children = true;
        subexpressions.push({to_split.subspan(offset, size)});
      }
      offset += size;
    }
    if (!has_children) {
      // There's nothing else to split.
      subexpressions.push({to_split});
      break;
    }
  }

  std::vector<absl::Span<const Token>> subexpressions_vector;
  subexpressions_vector.reserve(subexpressions.size());
  while (!subexpressions.empty()) {
    subexpressions_vector.push_back(subexpressions.top().slice);
    subexpressions.pop();
  }

  // Sort all subexpressions according to their position in the global
  // expression.
  std::sort(subexpressions_vector.begin(), subexpressions_vector.end(),
            [](const absl::Span<const Token>& lhs,
               const absl::Span<const Token>& rhs) {
              return lhs.data() < rhs.data();
            });
  return subexpressions_vector;
}

// A subexpression and its partial solution.
struct SubexpressionSolution {
  // The expression to use as the constraint for this subproblem.
  //
  // A slice into the global expression.
  absl::Span<const Token> expression;

  // The items selected in the solution.
  std::vector<int> selected;

  // The partial dual cost of the solution.
  //
  // This is (c . x + lambda . (Weight x)).  In other words, it's the
  // component of the dual which is a function of 'x', the selected items.
  double partial_dual_cost;

  // Whether this subexpression has a successfully-computed solution.
  bool success;
};

}  // namespace

struct ParallelSolver::Workspace {
  // Precomputed table of subexpression sizes for the current SelectionProblem's
  // expression constraint.
  std::vector<int> subexpression_size;

  // All work items, sorted by their start index.
  std::vector<SubexpressionSolution> subexpression_solutions;
};

ParallelSolver::ParallelSolver(int thread_count, SolverFactory solver_factory)
    : workspace_(new Workspace),
      thread_count_(thread_count),
      solver_factory_(std::move(solver_factory)) {}
ParallelSolver::~ParallelSolver() {}

void ParallelSolver::Init(const SelectionProblem& problem) {
  // Generate more subexpressions than threads to even out thread utilization.
  const int kSubexpressionsPerThread = 2;

  problem_ = problem;

  workspace_->subexpression_size.resize(problem_.expression.size());
  PrecomputeSubexpressionSize(problem_.expression,
                              absl::MakeSpan(workspace_->subexpression_size));

  int num_subexpressions = thread_count_ * kSubexpressionsPerThread;
  std::vector<absl::Span<const Token>> subexpressions = SplitIntoSubexpressions(
      problem_.expression, workspace_->subexpression_size, num_subexpressions);

  while (dual_solvers_.size() < subexpressions.size()) {
    dual_solvers_.push_back(solver_factory_());
  }

  // Initialize the solver for each subexpression.
  workspace_->subexpression_solutions.resize(subexpressions.size());
  base::BalancedParallelFor(thread_count_, subexpressions.size(), [&](int s) {
    DualSelectionSolver& dual_solver = *dual_solvers_[s];
    workspace_->subexpression_solutions[s].expression = subexpressions[s];
    SubexpressionSolution& work_item = workspace_->subexpression_solutions[s];
    SelectionProblem subproblem;
    subproblem.expression = work_item.expression;
    subproblem.capacity = problem_.capacity;
    subproblem.items = problem_.items;
    dual_solver.Init(subproblem);
  });
}

namespace {

// Recursively solves a selection problem, using precomputed solutions where
// possible.
//
// Performance here is not so critical because this will typically only solve
// the top of the expression tree, with large subtrees being solved in parallel
// to precompute SubexpressionSolutions.
bool SolveRecursive(
    const ItemSet& items, absl::Span<const Token> expression,
    absl::Span<const int> expression_size, absl::Span<const double> multipliers,
    absl::Span<const SubexpressionSolution> precomputed_solutions,
    std::vector<int>* solution_vars, double* partial_dual_cost) {
  if (expression.empty()) {
    return false;
  }
  const Token& token = expression.front();
  // If there is a precomputed solution, use it.
  if (token.IsAnd() || token.IsOr()) {
    const auto precomputed_solution = std::lower_bound(
        precomputed_solutions.begin(), precomputed_solutions.end(),
        expression.data(),
        [](const SubexpressionSolution& s, const Token* token_ptr) {
          return s.expression.begin() < token_ptr;
        });
    if (precomputed_solution->expression.data() == expression.data()) {
      solution_vars->insert(solution_vars->end(),
                            precomputed_solution->selected.begin(),
                            precomputed_solution->selected.end());
      *partial_dual_cost = precomputed_solution->partial_dual_cost;
      return precomputed_solution->success;
    }
  }
  // Otherwise, solve recursively.
  if (token.IsAnd()) {
    *partial_dual_cost = 0.0;
    for (int expression_index = 1; expression_index < expression.size() - 1;
         ++expression_index) {
      if (expression[expression_index].IsEnd()) break;
      int size = expression_size[expression_index];
      double subexpression_partial_dual_cost;
      if (!SolveRecursive(items, expression.subspan(expression_index, size),
                          expression_size.subspan(expression_index, size),
                          multipliers, precomputed_solutions, solution_vars,
                          &subexpression_partial_dual_cost)) {
        return false;
      }
      *partial_dual_cost += subexpression_partial_dual_cost;
      expression_index += size - 1;
    }
    return true;
  } else if (token.IsOr()) {
    bool any_success = false;
    *partial_dual_cost = 0.0;
    const int solution_vars_begin = solution_vars->size();
    for (int expression_index = 1; expression_index < expression.size() - 1;
         ++expression_index) {
      if (expression[expression_index].IsEnd()) break;
      int solution_vars_head = solution_vars->size();
      int size = expression_size[expression_index];
      double subexpression_partial_dual_cost;
      bool subexpression_success =
          SolveRecursive(items, expression.subspan(expression_index, size),
                         expression_size.subspan(expression_index, size),
                         multipliers, precomputed_solutions, solution_vars,
                         &subexpression_partial_dual_cost);
      if (subexpression_success &&
          (!any_success ||
           subexpression_partial_dual_cost < *partial_dual_cost)) {
        // Use this solution.
        *partial_dual_cost = subexpression_partial_dual_cost;
        std::copy(solution_vars->begin() + solution_vars_head,
                  solution_vars->end(),
                  solution_vars->begin() + solution_vars_begin);
        solution_vars->resize(solution_vars->size() -
                              (solution_vars_head - solution_vars_begin));
        any_success = true;
      } else {
        solution_vars->resize(solution_vars_head);
      }
      expression_index += size - 1;
    }
  } else {
    DCHECK(token.IsItem());
    int item = token.GetItem();
    *partial_dual_cost = items.GetCost(item);
    for (const auto& weight : items.GetWeights(item)) {
      *partial_dual_cost += multipliers[weight.index] * weight.value;
    }
    solution_vars->push_back(item);
  }
  return true;
}

}  // namespace

bool ParallelSolver::Solve(absl::Span<const double> multipliers,
                           std::vector<int>* selected) {
  DCHECK_EQ(problem_.capacity.size(), multipliers.size());
  selected->clear();
  // Solve problems corresponding to the disjoint subexpressions in parallel.
  base::BalancedParallelFor(
      thread_count_, workspace_->subexpression_solutions.size(), [&](int s) {
        DualSelectionSolver& dual_solver = *dual_solvers_[s];
        SubexpressionSolution& work_item =
            workspace_->subexpression_solutions[s];
        work_item.success = dual_solver.Solve(multipliers, &work_item.selected);

        // Compute the partial dual cost for the solution to this subexpression.
        //
        // TODO(puneetl):  In practice, the RelaxedSolver instance which this
        // typically uses for dual_solver would have already computed a
        // partial-dual-cost.  It may be worth optimizing to reuse that value.
        work_item.partial_dual_cost = 0.0;
        for (int item : work_item.selected) {
          work_item.partial_dual_cost += problem_.items->GetCost(item);
          for (const auto& weight : problem_.items->GetWeights(item)) {
            work_item.partial_dual_cost +=
                multipliers[weight.index] * weight.value;
          }
        }
      });

  // Solve the global problem, using the precomputed subexpression solutions.
  double partial_dual_cost = 0.0f;
  return SolveRecursive(*problem_.items, problem_.expression,
                        workspace_->subexpression_size, multipliers,
                        workspace_->subexpression_solutions, selected,
                        &partial_dual_cost);
}

}  // namespace selection
}  // namespace tiler
}  // namespace seurat
