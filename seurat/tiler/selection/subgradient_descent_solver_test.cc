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

#include "seurat/tiler/selection/subgradient_descent_solver.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/tiler/selection/bisection_solver.h"
#include "seurat/tiler/selection/relaxed_solver.h"
#include "seurat/tiler/selection/selection_test_util.h"
#include "seurat/tiler/selection/sequential_solver.h"

namespace seurat {
namespace tiler {
namespace selection {
namespace {

std::unique_ptr<SelectionSolver> MakeSolverToTest() {
  const int kThreadCount = 1;
  const int kBisectionIterations = 100;
  const int kSubgradientDescentIterations = 1000;
  std::vector<std::unique_ptr<SelectionSolver>> sequence;
  sequence.emplace_back(new SubgradientDescentSolver(
      kThreadCount, std::unique_ptr<DualSelectionSolver>(new RelaxedSolver),
      kSubgradientDescentIterations));
  sequence.emplace_back(new BisectionSolver(
      kThreadCount, std::unique_ptr<DualSelectionSolver>(new RelaxedSolver),
      kBisectionIterations));
  return std::unique_ptr<SelectionSolver>(
      new SequentialSolver(std::move(sequence)));
}

TEST(SubgradientDescentSolverTest, Test) {
  TestSolverFindsReasonableSolution(MakeSolverToTest().get());
}

TEST(SubgradientDescentSolverTest, Test_NoItems) {
  TestSolverFailsWhenNoItemsPresent(MakeSolverToTest().get());
}

TEST(SubgradientDescentSolverTest, Test_SingleItem) {
  TestSolverWithSingleItem(MakeSolverToTest().get());
}

}  // namespace
}  // namespace selection
}  // namespace tiler
}  // namespace seurat
