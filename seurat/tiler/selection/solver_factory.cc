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

#include "seurat/tiler/selection/solver_factory.h"

#include "seurat/tiler/selection/bisection_solver.h"
#include "seurat/tiler/selection/parallel_solver.h"
#include "seurat/tiler/selection/relaxed_solver.h"
#include "seurat/tiler/selection/selection_util.h"
#include "seurat/tiler/selection/sequential_solver.h"
#include "seurat/tiler/selection/subgradient_descent_solver.h"

namespace seurat {
namespace tiler {
namespace selection {

std::unique_ptr<SelectionSolver> CreateSelectionSolver(
    const SelectionSolverParameters& params) {
  auto make_dual_solver = [=]() {
    return std::unique_ptr<DualSelectionSolver>(
        new ParallelSolver(params.thread_count, []() {
          return std::unique_ptr<DualSelectionSolver>(new RelaxedSolver);
        }));
  };
  std::vector<std::unique_ptr<SelectionSolver>> sequence;
  // Start by bisecting to try to get an initial solution.
  sequence.emplace_back(new BisectionSolver(
      params.thread_count, make_dual_solver(), params.bisection_iterations,
      /* initialize multipliers = */ true));
  // Optimize with subgradient descent to get a high-quality (but not
  // necessarily feasible!) solution.
  sequence.emplace_back(
      new SubgradientDescentSolver(params.thread_count, make_dual_solver(),
                                   params.subgradient_descent_iterations));
  // Refine the result with bisection to get a feasible solution.
  sequence.emplace_back(
      new BisectionSolver(params.thread_count, make_dual_solver(),
                          params.bisection_iterations, false));
  return std::unique_ptr<SelectionSolver>(
      new SequentialSolver(std::move(sequence)));
}

}  // namespace selection
}  // namespace tiler
}  // namespace seurat
