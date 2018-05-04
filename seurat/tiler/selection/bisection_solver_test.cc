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

#include "seurat/tiler/selection/bisection_solver.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/tiler/selection/relaxed_solver.h"
#include "seurat/tiler/selection/selection_test_util.h"

namespace seurat {
namespace tiler {
namespace selection {
namespace {

std::unique_ptr<SelectionSolver> MakeSolverToTest() {
  const int kBisectionIterations = 100;
  const bool kInitializeMultipliers = true;
  std::unique_ptr<SelectionSolver> solver(new BisectionSolver(
      1, std::unique_ptr<DualSelectionSolver>(new RelaxedSolver),
      kBisectionIterations, kInitializeMultipliers));
  return solver;
}

TEST(BisectionSolverTest, TestBisection) {
  TestSolverFindsReasonableSolution(MakeSolverToTest().get());
}

TEST(BisectionSolverTest, TestBisection_NoItems) {
  TestSolverFailsWhenNoItemsPresent(MakeSolverToTest().get());
}

TEST(BisectionSolverTest, TestBisection_SingleItem) {
  TestSolverWithSingleItem(MakeSolverToTest().get());
}

}  // namespace
}  // namespace selection
}  // namespace tiler
}  // namespace seurat
