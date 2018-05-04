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

#include "seurat/tiler/selection/selection_test_util.h"

#include <array>
#include <random>
#include <vector>

#include "gtest/gtest.h"
#include "absl/types/span.h"
#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_util.h"

namespace seurat {
namespace tiler {
namespace selection {

namespace {

bool ItemInSelection(absl::Span<const int> selected_items, int item) {
  return std::find(selected_items.begin(), selected_items.end(), item) !=
         selected_items.end();
}

}  // namespace

void TestSolverFindsReasonableSolution(SelectionSolver *solver) {
  // Solve a problem consisting of the following expression:
  //  (a || c || f) && (b || ((d || g) && (e || (h && i && j))))
  //
  // Equivalently:
  //  AND(OR(a, c, f), OR(b, AND(OR(d, g), OR(e, AND(h, i, j)))))

  std::mt19937 random;
  std::uniform_real_distribution<float> cost(0.1f, 10.0f);
  std::uniform_real_distribution<float> weight0(1.0f, 10.0f);
  std::uniform_real_distribution<float> weight1(20.0f, 50.0f);

  ItemSet items(2);

  // Helper lambda to create an item with the specified |cost| and a random
  // weight.
  auto append_item_with_random_weight = [&](float cost) -> int {
    return items.AppendItem(cost,
                            {{{0, weight0(random)}, {1, weight1(random)}}});
  };

  // Make 'a', 'b', and 'i' have such high cost that they should certainly not
  // be in the final solution.

  const float kExpensive = 1.0e5f;

  int a = append_item_with_random_weight(kExpensive);
  int b = append_item_with_random_weight(kExpensive);
  int c = append_item_with_random_weight(cost(random));
  int d = append_item_with_random_weight(cost(random));
  int e = append_item_with_random_weight(cost(random));
  int f = append_item_with_random_weight(cost(random));
  int g = append_item_with_random_weight(cost(random));
  int h = append_item_with_random_weight(cost(random));
  int i = append_item_with_random_weight(kExpensive);
  int j = append_item_with_random_weight(cost(random));

  std::vector<Token> expression = {
      Token::And(),   Token::Or(),    Token::Item(a), Token::Item(c),
      Token::Item(f), Token::End(),   Token::Or(),    Token::Item(b),
      Token::And(),   Token::Or(),    Token::Item(d), Token::Item(g),
      Token::End(),   Token::Or(),    Token::Item(e), Token::And(),
      Token::Item(h), Token::Item(i), Token::Item(j), Token::End(),
      Token::End(),   Token::End(),   Token::End(),   Token::End()};

  // Note that the capacity is large enough that it will not actively constrain
  // the solution.
  std::vector<double> capacity = {{10.0 * 10, 50.0 * 10}};

  SelectionProblem problem;
  problem.items = &items;
  problem.expression = expression;
  problem.capacity = capacity;

  std::vector<int> solution_items;
  std::vector<double> solution_weight_multipliers(2);

  solver->Init(problem);
  EXPECT_TRUE(solver->Solve(absl::MakeSpan(solution_weight_multipliers),
                            &solution_items));

  // Verify that whatever solution was generated is feasible.

  int a_s = ItemInSelection(solution_items, a) ? 1 : 0;
  int b_s = ItemInSelection(solution_items, b) ? 1 : 0;
  int c_s = ItemInSelection(solution_items, c) ? 1 : 0;
  int d_s = ItemInSelection(solution_items, d) ? 1 : 0;
  int e_s = ItemInSelection(solution_items, e) ? 1 : 0;
  int f_s = ItemInSelection(solution_items, f) ? 1 : 0;
  int g_s = ItemInSelection(solution_items, g) ? 1 : 0;
  int h_s = ItemInSelection(solution_items, h) ? 1 : 0;
  int i_s = ItemInSelection(solution_items, i) ? 1 : 0;
  int j_s = ItemInSelection(solution_items, j) ? 1 : 0;

  EXPECT_TRUE((a_s || c_s || f_s) &&
              (b_s || ((d_s || g_s) && (e_s || (h_s && i_s && j_s)))));
  std::vector<double> solution_weight(2);
  ComputeTotalWeight(items, solution_items, absl::MakeSpan(solution_weight));
  EXPECT_TRUE(IsFeasibleWeight(problem, solution_weight));
}

void TestSolverFailsWhenNoItemsPresent(SelectionSolver *solver) {
  // Test with various expressions, each of which has no actual items.
  ItemSet items(1);
  std::vector<double> capacity = {1.0};
  std::vector<int> selected;
  std::vector<double> multipliers = {1.0};
  {
    // And()
    std::vector<Token> expression = {Token::And(), Token::End()};
    solver->Init({&items, expression, capacity});
    EXPECT_FALSE(solver->Solve(absl::MakeSpan(multipliers), &selected));
  }

  {
    // Or()
    std::vector<Token> expression = {Token::Or(), Token::End()};
    solver->Init({&items, expression, capacity});
    EXPECT_FALSE(solver->Solve(absl::MakeSpan(multipliers), &selected));
  }

  {
    // And(And())
    std::vector<Token> expression = {Token::And(),  //
                                     Token::And(),  //
                                     Token::End(),  //
                                     Token::End()};
    solver->Init({&items, expression, capacity});
    EXPECT_FALSE(solver->Solve(absl::MakeSpan(multipliers), &selected));
  }

  {
    // And(And(), Or(), And(And(Or())))
    std::vector<Token> expression = {
        Token::And(),  //
        //
        Token::And(),  //
        Token::End(),  //
        //
        Token::Or(),   //
        Token::End(),  //
        //
        Token::And(),  //
        Token::And(),  //
        Token::Or(),   //
        Token::End(),  //
        Token::End(),  //
        Token::End(),  //
        //
        Token::End()  //
    };
    solver->Init({&items, expression, capacity});
    EXPECT_FALSE(solver->Solve(absl::MakeSpan(multipliers), &selected));
  }
}

void TestSolverWithSingleItem(SelectionSolver *solver) {
  // Construct a problem with a single item in a convoluted expression:
  // Or(Or(And(0)), And())
  //
  // All solvers should be able to find that the solution is the item set {0}.
  ItemSet items(3);
  std::vector<double> capacity = {2.0, 10.0, 4.0};
  const int item = items.AppendItem(
      4.0f, std::array<ItemSet::Weight, 2>{
                {ItemSet::Weight{0, 2.0}, ItemSet::Weight{2, 4.0}}});

  std::vector<double> multipliers = {2.0, 100.0, 4.0};
  std::vector<int> selected;

  std::vector<Token> expression = {
      Token::Or(),   //
      Token::Or(),   //
      Token::And(),  //
      Token::Item(item),
      Token::End(),  //
      Token::End(),  //
      Token::And(),  //
      Token::End(),  //
      Token::End(),  //
  };
  solver->Init({&items, expression, capacity});
  EXPECT_TRUE(solver->Solve(absl::MakeSpan(multipliers), &selected));
  EXPECT_EQ(1, selected.size());
  EXPECT_EQ(item, selected[0]);
}

}  // namespace selection
}  // namespace tiler
}  // namespace seurat
