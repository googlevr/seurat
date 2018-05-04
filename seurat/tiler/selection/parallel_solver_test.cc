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
  // It's infeasible to test the ParallelSolver in isolation.
  //
  // Instead, test it in the context of a minimal, complete solver which only
  // performs bisection.
  const int kBisectionIterations = 100;
  return std::unique_ptr<SelectionSolver>(new BisectionSolver(
      thread_count,
      std::unique_ptr<DualSelectionSolver>(new ParallelSolver(
          thread_count,
          []() {
            return std::unique_ptr<DualSelectionSolver>(new RelaxedSolver);
          })),
      kBisectionIterations));
}

TEST(ParallelSolverTest, TestSingleThreadNoItems) {
  TestSolverFailsWhenNoItemsPresent(MakeSolverToTest(1).get());
}

TEST(ParallelSolverTest, TestSingleThreadSingleItem) {
  TestSolverWithSingleItem(MakeSolverToTest(1).get());
}

TEST(ParallelSolverTest, TestSingleThreadMultipleItems) {
  TestSolverFindsReasonableSolution(MakeSolverToTest(1).get());
}

TEST(ParallelSolverTest, TestMultipleThreadsNoItems) {
  TestSolverFailsWhenNoItemsPresent(MakeSolverToTest(3).get());
}

TEST(ParallelSolverTest, TestMultipleThreadsSingleItem) {
  TestSolverWithSingleItem(MakeSolverToTest(3).get());
}

TEST(ParallelSolverTest, TestMultipleThreadsMultipleItems) {
  TestSolverFindsReasonableSolution(MakeSolverToTest(3).get());
}

}  // namespace
}  // namespace selection
}  // namespace tiler
}  // namespace seurat
