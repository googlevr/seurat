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

#include <array>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/tiler/selection/selection_problem.h"

namespace seurat {
namespace tiler {
namespace selection {
namespace {

TEST(SelectionUtilTest, TestIsFeasibleWeight) {
  std::vector<double> capacity = {2.0, 3.0, 1.0};
  EXPECT_TRUE(IsFeasibleWeight(SelectionProblem{nullptr, {}, capacity},
                               std::vector<double>{0.0, 2.9, 0.5}));
  EXPECT_FALSE(IsFeasibleWeight(SelectionProblem{nullptr, {}, capacity},
                                std::vector<double>{0.0, 3.1, 0.5}));
}

TEST(SelectionUtilTest, TestComputeTotalCost) {
  ItemSet items(0);
  int item3 = items.AppendItem(3.0, std::array<ItemSet::Weight, 0>{{}});
  int item2 = items.AppendItem(2.0, std::array<ItemSet::Weight, 0>{{}});
  int item9 = items.AppendItem(9.0, std::array<ItemSet::Weight, 0>{{}});

  EXPECT_EQ(3.0, ComputeTotalCost(items, std::array<int, 1>{{item3}}));
  EXPECT_EQ(2.0, ComputeTotalCost(items, std::array<int, 1>{{item2}}));
  EXPECT_EQ(9.0, ComputeTotalCost(items, std::array<int, 1>{{item9}}));
  EXPECT_NEAR(
      3.0 + 2.0 + 9.0,
      ComputeTotalCost(items, std::array<int, 3>{{item3, item2, item9}}), 1e-6);
  EXPECT_NEAR(2.0 + 9.0,
              ComputeTotalCost(items, std::array<int, 2>{{item2, item9}}),
              1e-6);
}

TEST(SelectionUtilTest, TestPrintExpressionToString) {
  // Build an expression for testing:
  //   AND(1, 2, 3, OR(4))
  std::vector<Token> expression = {
      Token::And(), Token::Item(1), Token::Item(2), Token::Item(3),
      Token::Or(),  Token::Item(4), Token::End(),   Token::End()};
  std::string str = PrintExpressionToString(expression);

  EXPECT_EQ("AND( 1 2 3 OR( 4 ) ) ", str);
}

TEST(SelectionUtilTest, TestPrecomputeSubexpressionSize) {
  // The expression:
  //   AND(1, 2, 3, OR(4))
  // should have subexpression sizes:
  //   8 1 1 1 3 1 3 8
  std::vector<Token> expression = {
      Token::And(), Token::Item(1), Token::Item(2), Token::Item(3),
      Token::Or(),  Token::Item(4), Token::End(),   Token::End()};

  std::vector<int> sizes(expression.size());
  PrecomputeSubexpressionSize(expression, absl::MakeSpan(sizes));

  std::vector<int> expected_sizes = {8, 1, 1, 1, 3, 1, 3, 8};
  EXPECT_EQ(expected_sizes, sizes);
}

TEST(SelectionUtilTest, TestValidateSelection) {
  {
    // AND(1, 2, 3, OR(4))
    std::vector<Token> expression = {
        Token::And(), Token::Item(1), Token::Item(2), Token::Item(3),
        Token::Or(),  Token::Item(4), Token::End(),   Token::End()};

    EXPECT_FALSE(ValidateSelection(expression, {}));
    EXPECT_FALSE(ValidateSelection(expression, std::vector<int>{1}));
    EXPECT_FALSE(ValidateSelection(expression, std::vector<int>{4}));
    EXPECT_FALSE(ValidateSelection(expression, std::vector<int>{1, 2, 3}));
    EXPECT_TRUE(ValidateSelection(expression, std::vector<int>{1, 2, 3, 4}));
  }

  {
    // OR( AND(), 1, OR())
    std::vector<Token> expression = {Token::Or(), Token::And(), Token::End(),
                                     Token::Item(1), Token::End()};

    EXPECT_FALSE(ValidateSelection(expression, {}));
    EXPECT_TRUE(ValidateSelection(expression, std::vector<int>{1}));
  }

  {
    // OR(AND(OR(), 1), 2, AND(3, OR()), AND(AND(), 4), AND(5, AND()))
    std::vector<Token> expression = {
        Token::Or(),    Token::And(), Token::Or(),    Token::End(),
        Token::Item(1), Token::End(), Token::Item(2), Token::And(),
        Token::Item(3), Token::Or(),  Token::End(),   Token::End(),
        Token::And(),   Token::And(), Token::End(),   Token::Item(4),
        Token::End(),   Token::And(), Token::Item(5), Token::And(),
        Token::End(),   Token::End(), Token::End()};

    // Selecting '2' is valid, but all others are not, since they would also
    // imply the use of an empty-subexpression.
    EXPECT_FALSE(ValidateSelection(expression, {}));
    EXPECT_FALSE(ValidateSelection(expression, std::vector<int>{1}));
    EXPECT_TRUE(ValidateSelection(expression, std::vector<int>{2}));
    EXPECT_FALSE(ValidateSelection(expression, std::vector<int>{3}));
    EXPECT_FALSE(ValidateSelection(expression, std::vector<int>{4}));
    EXPECT_FALSE(ValidateSelection(expression, std::vector<int>{5}));
  }
}

}  // namespace
}  // namespace selection
}  // namespace tiler
}  // namespace seurat
