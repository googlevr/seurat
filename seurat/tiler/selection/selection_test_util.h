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

#ifndef VR_SEURAT_TILER_SELECTION_SELECTION_TEST_UTIL_H_
#define VR_SEURAT_TILER_SELECTION_SELECTION_TEST_UTIL_H_

#include "seurat/tiler/selection/selection_solver.h"

namespace seurat {
namespace tiler {
namespace selection {

// Runs the solver on a small example problem and verifies that the generated
// solution is feasible & not obviously suboptimal.
void TestSolverFindsReasonableSolution(SelectionSolver* solver);

// Runs the solver on various problems consisting of expression constraints
// which do not have any items.
void TestSolverFailsWhenNoItemsPresent(SelectionSolver* solver);

// Runs the solver on a problem with only 1 relevant item.
void TestSolverWithSingleItem(SelectionSolver* solver);

}  // namespace selection
}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SELECTION_SELECTION_TEST_UTIL_H_
