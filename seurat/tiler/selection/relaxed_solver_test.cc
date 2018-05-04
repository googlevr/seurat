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

#include <array>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_test_util.h"
#include "seurat/tiler/selection/selection_util.h"

namespace seurat {
namespace tiler {
namespace selection {
namespace {

TEST(RelaxedSolverTest, TestNoItems) {
  ItemSet items(1);
  std::vector<double> capacity = {1.0};
  RelaxedSolver solver;

  std::vector<int> result;
  std::vector<double> multipliers = {1.0};

  {
    std::vector<Token> expression = {Token::And(), Token::End()};
    solver.Init({&items, expression, capacity});
    EXPECT_FALSE(solver.Solve(multipliers, &result));
  }

  {
    std::vector<Token> expression = {Token::Or(), Token::End()};
    solver.Init({&items, expression, capacity});
    EXPECT_FALSE(solver.Solve(multipliers, &result));
  }

  {
    std::vector<Token> expression = {Token::And(),  //
                                     Token::And(),  //
                                     Token::End(),  //
                                     Token::End()};
    solver.Init({&items, expression, capacity});
    EXPECT_FALSE(solver.Solve(multipliers, &result));
  }

  {
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
    solver.Init({&items, expression, capacity});
    EXPECT_FALSE(solver.Solve(multipliers, &result));
  }
}

TEST(RelaxedSolverTest, TestSingleItem) {
  ItemSet items(3);
  std::vector<double> capacity = {2.0, 10.0, 4.0};
  const int item = items.AppendItem(
      4.0f, std::array<ItemSet::Weight, 2>{
                {ItemSet::Weight{0, 2.0}, ItemSet::Weight{2, 4.0}}});

  RelaxedSolver solver;

  std::vector<double> multipliers = {2.0, 100.0, 4.0};
  std::vector<int> result;

  // Or(Or(And(i)), And())
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
  solver.Init({&items, expression, capacity});
  EXPECT_TRUE(solver.Solve(multipliers, &result));
  EXPECT_EQ(1, result.size());
  EXPECT_EQ(item, result[0]);
}

TEST(RelaxedSolverTest, TestSingleItem_2) {
  ItemSet items(3);
  std::vector<double> capacity = {2.0, 10.0, 4.0};
  const int item = items.AppendItem(
      4.0f, std::array<ItemSet::Weight, 2>{
                {ItemSet::Weight{0, 2.0}, ItemSet::Weight{2, 4.0}}});

  RelaxedSolver solver;

  std::vector<double> multipliers = {2.0, 100.0, 4.0};
  std::vector<int> result;

  // Or(Or(And()), And(i)))
  std::vector<Token> expression = {
      Token::Or(),   //
      Token::Or(),   //
      Token::And(),  //
      Token::End(),  //
      Token::End(),  //
      Token::And(),  //
      Token::Item(item),
      Token::End(),  //
      Token::End(),  //
  };
  solver.Init({&items, expression, capacity});
  EXPECT_TRUE(solver.Solve(multipliers, &result));
  EXPECT_EQ(1, result.size());
  EXPECT_EQ(item, result[0]);
}

TEST(RelaxedSolverTest, TestSingleItem_Fail) {
  ItemSet items(3);
  std::vector<double> capacity = {2.0, 10.0, 4.0};
  const int item0 = items.AppendItem(
      4.0f, std::array<ItemSet::Weight, 2>{
                {ItemSet::Weight{0, 2.0}, ItemSet::Weight{2, 4.0}}});
  const int item1 = items.AppendItem(
      4.0f, std::array<ItemSet::Weight, 2>{
                {ItemSet::Weight{0, 2.0}, ItemSet::Weight{2, 4.0}}});

  RelaxedSolver solver;

  std::vector<double> multipliers = {2.0, 100.0, 4.0};
  std::vector<int> result;

  // And(Or(Or(), And(Or(), 0, Or())), 1)
  std::vector<Token> expression = {
      Token::And(),        //
      Token::Or(),         //
      Token::Or(),         //
      Token::End(),        //
      Token::And(),        //
      Token::Or(),         //
      Token::End(),        //
      Token::Item(item0),  //
      Token::Or(),         //
      Token::End(),        //
      Token::End(),        //
      Token::End(),        //
      Token::Item(item1),  //
      Token::End()         //
  };
  solver.Init({&items, expression, capacity});
  EXPECT_FALSE(solver.Solve(multipliers, &result));
}

}  // namespace
}  // namespace selection
}  // namespace tiler
}  // namespace seurat
