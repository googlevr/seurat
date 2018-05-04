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

#include "seurat/tiler/selection/selection_util.h"

#include "ion/base/logging.h"
#include "seurat/base/parallel.h"

namespace seurat {
namespace tiler {
namespace selection {

bool IsFeasibleWeight(const SelectionProblem& problem,
                      const absl::Span<const double> total_weight) {
  DCHECK_EQ(total_weight.size(), problem.capacity.size());
  for (int d = 0; d < total_weight.size(); ++d) {
    if (total_weight[d] > problem.capacity[d]) {
      return false;
    }
  }
  return true;
}

double ComputeTotalCost(const ItemSet& item_set, absl::Span<const int> items) {
  double total_cost = 0.0;
  for (int i : items) {
    total_cost += item_set.GetCost(i);
  }
  return total_cost;
}

void ComputeTotalWeight(const ItemSet& item_set, absl::Span<const int> items,
                        absl::Span<double> weight) {
  std::fill(weight.begin(), weight.end(), 0.0);
  for (int i : items) {
    for (const auto& w : item_set.GetWeights(i)) {
      weight[w.index] += w.value;
    }
  }
}

double ComputeDualCost(const SelectionProblem& problem,
                       absl::Span<const double> weight_multipliers,
                       absl::Span<const double> weight, double cost) {
  double dual_cost = cost;
  for (int w = 0; w < weight.size(); ++w) {
    dual_cost += weight[w] * weight_multipliers[w];
    dual_cost -= weight_multipliers[w] * problem.capacity[w];
  }
  return dual_cost;
}

std::string PrintExpressionToString(absl::Span<const Token> expression) {
  std::string str;
  for (const auto& token : expression) {
    if (token.IsAnd()) {
      str += "AND( ";
    } else if (token.IsOr()) {
      str += "OR( ";
    } else if (token.IsEnd()) {
      str += ") ";
    } else {
      str += std::to_string(token.GetItem());
      str += " ";
    }
  }
  return str;
}

void PrecomputeSubexpressionSize(absl::Span<const Token> expression,
                                 absl::Span<int> subexpression_size) {
  CHECK_EQ(expression.size(), subexpression_size.size());
  std::vector<int> stack;
  for (int i = 0; i < expression.size(); ++i) {
    Token token = expression[i];
    if (token.IsAnd() || token.IsOr()) {
      stack.push_back(i);
    } else if (token.IsEnd()) {
      int start_index = stack.back();
      stack.pop_back();
      subexpression_size[i] = (i - start_index) + 1;
      subexpression_size[start_index] = (i - start_index) + 1;
    } else {
      CHECK(token.IsItem());
      subexpression_size[i] = 1;
    }
  }
  CHECK(stack.empty()) << "Unbalanced expression";
}

namespace {

bool ValidateSelectionRecursive(absl::Span<const Token> expression,
                                const std::vector<int>& sorted_selection,
                                int* tokens_consumed) {
  *tokens_consumed = 1;
  const Token& token = expression.front();
  if (token.IsAnd()) {
    int offset = 1;
    bool success = true;
    while (!expression[offset].IsEnd()) {
      int subexpression_size;
      if (!ValidateSelectionRecursive(
              expression.subspan(offset, expression.size() - offset),
              sorted_selection, &subexpression_size)) {
        success = false;
      }
      *tokens_consumed += subexpression_size;
      offset += subexpression_size;
    }
    // Consumed the End() token.
    *tokens_consumed += 1;
    if (offset == 1) {
      // The empty And() clause always evaluates to false.
      success = false;
    }
    return success;
  } else if (token.IsOr()) {
    int offset = 1;
    // Empty Or() clause always evaluates to false.
    bool success = false;
    while (!expression[offset].IsEnd()) {
      int subexpression_size;
      if (ValidateSelectionRecursive(
              expression.subspan(offset, expression.size() - offset),
              sorted_selection, &subexpression_size)) {
        success = true;
      }
      *tokens_consumed += subexpression_size;
      offset += subexpression_size;
    }
    // Consume the End() token.
    *tokens_consumed += 1;
    return success;
  } else {
    CHECK(token.IsItem());
    int item = token.GetItem();
    *tokens_consumed = 1;
    return std::binary_search(sorted_selection.begin(), sorted_selection.end(),
                              item);
  }
}

}  // namespace

bool ValidateSelection(absl::Span<const Token> expression,
                       absl::Span<const int> selected_items) {
  std::vector<int> sorted_items(selected_items.begin(), selected_items.end());
  std::sort(sorted_items.begin(), sorted_items.end());
  int unused;
  return ValidateSelectionRecursive(expression, sorted_items, &unused);
}

}  // namespace selection
}  // namespace tiler
}  // namespace seurat
