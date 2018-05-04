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

#include "seurat/tiler/selection/cost_calculator.h"

#include <array>
#include <random>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_util.h"

namespace seurat {
namespace tiler {
namespace selection {
namespace {

TEST(CostCalculatorTest, TestComputeSolutionCost) {
  std::mt19937 random;
  std::uniform_real_distribution<float> cost(0.1f, 10.0f);
  std::uniform_real_distribution<float> weight0(1.0f, 10.0f);
  std::uniform_real_distribution<float> weight1(20.0f, 50.0f);

  ItemSet items(2);

  const int kNumItems = 10;
  for (int i = 0; i < kNumItems; ++i) {
    items.AppendItem(cost(random),
                     {{{0, weight0(random)}, {1, weight1(random)}}});
  }

  std::vector<double> capacity = {{10.0 * 10, 50.0 * 10}};

  SelectionProblem problem;
  problem.items = &items;
  problem.capacity = capacity;

  // Test with various thread counts.
  for (int thread_count :
       std::array<int, 11>{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20}}) {
    CostCalculator cost_calculator(thread_count);

    std::vector<int> selected_items = {0, 1, 3, 6, 9};
    double primal_cost;
    double dual_cost;
    std::vector<double> multipliers = {30.0f, 42.0f};
    std::vector<double> total_weight(2);
    cost_calculator.ComputeSolutionCost(problem, multipliers, selected_items,
                                        &primal_cost, &dual_cost,
                                        absl::MakeSpan(total_weight));

    EXPECT_EQ(primal_cost, ComputeTotalCost(items, selected_items));
    EXPECT_EQ(dual_cost,
              ComputeDualCost(problem, multipliers, total_weight, primal_cost));
    std::vector<double> total_weight_2(2);
    ComputeTotalWeight(items, selected_items, absl::MakeSpan(total_weight_2));
    EXPECT_EQ(total_weight_2, total_weight);
  }
}

}  // namespace
}  // namespace selection
}  // namespace tiler
}  // namespace seurat
