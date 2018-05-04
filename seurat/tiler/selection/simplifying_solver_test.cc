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

#include <array>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/tiler/selection/bisection_solver.h"
#include "seurat/tiler/selection/relaxed_solver.h"
#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_test_util.h"
#include "seurat/tiler/selection/selection_util.h"

namespace seurat {
namespace tiler {
namespace selection {
namespace {

std::unique_ptr<SelectionSolver> MakeSolverToTest(int thread_count) {
  // It's infeasible to test the SimplifyingSolver in isolation.
  //
  // Instead, test it in the context of a minimal, complete solver which only
  // performs bisection.
  const int kBisectionIterations = 100;
  return std::unique_ptr<SelectionSolver>(new SimplifyingSolver(
      std::unique_ptr<SelectionSolver>(new BisectionSolver(
          thread_count, std::unique_ptr<DualSelectionSolver>(new RelaxedSolver),
          kBisectionIterations))));
}

TEST(SimplifyingSolverTest, TestSimpleExpressions) {
  ItemSet items(0);

  int item1 = items.AppendItem(10.0f, {});
  int item2 = items.AppendItem(10.0f, {});
  int item3 = items.AppendItem(1.0f, {});
  int item4 = items.AppendItem(1.0f, {});
  int item5 = items.AppendItem(1.0f, {});
  // Test with the expression constraint:
  //   Or(1, 2, And(3, 4, 5))
  // The result should be {3, 4, 5} with total cost 3.0.
  std::vector<Token> expression;
  expression.push_back(Token::Or());
  expression.push_back(Token::Item(item1));
  expression.push_back(Token::Item(item2));
  expression.push_back(Token::And());
  expression.push_back(Token::Item(item3));
  expression.push_back(Token::Item(item4));
  expression.push_back(Token::Item(item5));
  expression.push_back(Token::End());
  expression.push_back(Token::End());

  std::unique_ptr<SelectionSolver> solver = MakeSolverToTest(3);
  solver->Init({&items, expression, {}});
  std::vector<double> multipliers;
  std::vector<int> solution;
  EXPECT_TRUE(solver->Solve(absl::MakeSpan(multipliers), &solution));
  EXPECT_EQ(3, solution.size());
  EXPECT_THAT(solution, ::testing::UnorderedElementsAreArray(
                            std::array<int, 3>{{item3, item4, item5}}));
}

TEST(SimplifyingSolverTest, TestEmptySubexpressions) {
  ItemSet items(0);

  int item1 = items.AppendItem(10.0f, {});
  int item2 = items.AppendItem(10.0f, {});
  int item3 = items.AppendItem(1.0f, {});
  int item4 = items.AppendItem(1.0f, {});
  int item5 = items.AppendItem(1.0f, {});

  // Test with the expression constraint:
  //   And(1, Or(2, Or()), And(3, 4, Or (5, And())))
  // The result should be {1, 2, 3, 4, 5}.
  std::vector<Token> expression;
  expression.push_back(Token::And());
  expression.push_back(Token::Item(item1));
  expression.push_back(Token::Or());
  expression.push_back(Token::Item(item2));
  expression.push_back(Token::Or());
  expression.push_back(Token::End());
  expression.push_back(Token::End());
  expression.push_back(Token::And());
  expression.push_back(Token::Item(item3));
  expression.push_back(Token::Item(item4));
  expression.push_back(Token::Or());
  expression.push_back(Token::Item(item5));
  expression.push_back(Token::And());
  expression.push_back(Token::End());
  expression.push_back(Token::End());
  expression.push_back(Token::End());
  expression.push_back(Token::End());

  std::unique_ptr<SelectionSolver> solver = MakeSolverToTest(3);
  solver->Init({&items, expression, {}});
  std::vector<double> multipliers;
  std::vector<int> solution;
  solver->Solve(absl::MakeSpan(multipliers), &solution);
  EXPECT_THAT(solution, ::testing::UnorderedElementsAreArray(std::array<int, 5>{
                            {item1, item2, item3, item4, item5}}));
}

TEST(SimplifyingSolverTest, TestSingleThreadNoItems) {
  TestSolverFailsWhenNoItemsPresent(MakeSolverToTest(1).get());
}

TEST(SimplifyingSolverTest, TestSingleThreadSingleItem) {
  TestSolverWithSingleItem(MakeSolverToTest(1).get());
}

TEST(SimplifyingSolverTest, TestSingleThreadMultipleItems) {
  TestSolverFindsReasonableSolution(MakeSolverToTest(1).get());
}

TEST(SimplifyingSolverTest, TestMultipleThreadsNoItems) {
  TestSolverFailsWhenNoItemsPresent(MakeSolverToTest(3).get());
}

TEST(SimplifyingSolverTest, TestMultipleThreadsSingleItem) {
  TestSolverWithSingleItem(MakeSolverToTest(3).get());
}

TEST(SimplifyingSolverTest, TestMultipleThreadsMultipleItems) {
  TestSolverFindsReasonableSolution(MakeSolverToTest(3).get());
}

}  // namespace
}  // namespace selection
}  // namespace tiler
}  // namespace seurat
