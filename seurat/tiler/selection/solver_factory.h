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

#ifndef VR_SEURAT_TILER_SELECTION_SOLVER_FACTORY_H_
#define VR_SEURAT_TILER_SELECTION_SOLVER_FACTORY_H_

#include <memory>

#include "seurat/tiler/selection/selection_problem.h"
#include "seurat/tiler/selection/selection_solver.h"

namespace seurat {
namespace tiler {
namespace selection {

// Parameters to tune the SelectionSolver to build.
struct SelectionSolverParameters {
  int thread_count = 1;
  int bisection_iterations = 100;
  int subgradient_descent_iterations = 400;
};

// Builds a SelectionSolver.
std::unique_ptr<SelectionSolver> CreateSelectionSolver(
    const SelectionSolverParameters& params);

}  // namespace selection
}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SELECTION_SOLVER_FACTORY_H_
