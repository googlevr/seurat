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

#include "seurat/tiler/selection/selection_problem.h"

#include <array>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace seurat {
namespace tiler {
namespace selection {
namespace {

TEST(SelectionProblemTest, TestItemSet) {
  ItemSet items(3);
  EXPECT_EQ(3, items.GetNumWeights());
  int new_item = items.AppendItem(
      2.0f, std::array<ItemSet::Weight, 2>{{ItemSet::Weight{0, 1.0f},  //
                                            ItemSet::Weight{2, 0.5f}}});
  EXPECT_EQ(2.0f, items.GetCost(new_item));

  EXPECT_EQ(2, items.GetWeights(new_item).size());

  EXPECT_EQ(0, items.GetWeights(new_item)[0].index);
  EXPECT_EQ(1.0f, items.GetWeights(new_item)[0].value);
  EXPECT_EQ(2, items.GetWeights(new_item)[1].index);
  EXPECT_EQ(0.5f, items.GetWeights(new_item)[1].value);
}

TEST(SelectionProblemTest, TestToken) {
  EXPECT_TRUE(Token::And().IsAnd());
  EXPECT_TRUE(Token::Or().IsOr());
  EXPECT_TRUE(Token::End().IsEnd());
  EXPECT_TRUE(Token::Item(3).IsItem());

  Token token_and = Token::And();
  Token invalid;
  EXPECT_TRUE(token_and.IsAnd());
  EXPECT_FALSE(token_and.IsOr());
  EXPECT_FALSE(token_and.IsItem());
  EXPECT_FALSE(token_and.IsEnd());
  EXPECT_TRUE(token_and == Token::And());
  EXPECT_FALSE(token_and != Token::And());
  EXPECT_FALSE(invalid == token_and);

  const int kItem = 3;
  Token token_item = Token::Item(kItem);
  EXPECT_FALSE(token_item.IsAnd());
  EXPECT_FALSE(token_item.IsOr());
  EXPECT_FALSE(token_item.IsEnd());
  EXPECT_TRUE(token_item.IsItem());
  EXPECT_EQ(kItem, token_item.GetItem());
  EXPECT_FALSE(token_item == token_and);
  EXPECT_TRUE(token_item != token_and);
}

TEST(SelectionProblemTest, TestEmptyWeights) {
  // Add two items, where the last item has no weights.
  //
  // This verifies that GetWeights() will not try to access past the end of the
  // ItemSet's internal vector of weights for the last item.
  ItemSet items(3);

  int item0 = items.AppendItem(2.0f, {});

  EXPECT_EQ(0, items.GetWeights(item0).size());

  int item1 = items.AppendItem(
      2.0f, std::array<ItemSet::Weight, 2>{{ItemSet::Weight{0, 1.0f},  //
                                            ItemSet::Weight{2, 0.5f}}});
  EXPECT_EQ(2, items.GetWeights(item1).size());
}

}  // namespace
}  // namespace selection
}  // namespace tiler
}  // namespace seurat
