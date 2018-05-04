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

#ifndef VR_SEURAT_TILER_SELECTION_SELECTION_SOLVER_H_
#define VR_SEURAT_TILER_SELECTION_SELECTION_SOLVER_H_

#include <vector>

#include "absl/types/span.h"
#include "seurat/tiler/selection/selection_problem.h"

namespace seurat {
namespace tiler {
namespace selection {

// An interface for solving SelectionProblems.
//
// Implementations are not thread safe unless specified otherwise.
class SelectionSolver {
 public:
  virtual ~SelectionSolver() = default;

  // Sets the current problem to solve, possibly precomputing values required
  // for the subsequent Solve().
  //
  // A single Init() may be followed by multiple Solve()'s to iteratively
  // improve upon a solution.
  virtual void Init(const SelectionProblem& problem) = 0;

  // Solves the current problem, possibly making use of an existing
  // |multipliers|.
  //
  // Returns the items which were selected as well as the |multipliers|
  // corresponding to this solution.
  //
  // Returns false upon failure, e.g. due to numerical stability problems.
  virtual bool Solve(absl::Span<double> multipliers,
                     std::vector<int>* selected) = 0;
};

// Evaluates the dual of a SelectionProblem.
//
// See selection_problem.h for details.
class DualSelectionSolver {
 public:
  virtual ~DualSelectionSolver() = default;

  // Sets the current problem to solve, possibly precomputing values required
  // for the subsequent Solve().
  //
  // A single Init() may be followed by multiple Solve()'s to iteratively
  // improve upon a solution.
  virtual void Init(const SelectionProblem& problem) = 0;

  // Evaluates the dual of the current problem at the given |multipliers|.
  //
  // Returns false upon failure, e.g. due to numerical stability problems.
  virtual bool Solve(absl::Span<const double> multipliers,
                     std::vector<int>* selected) = 0;
};

}  // namespace selection
}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SELECTION_SELECTION_SOLVER_H_
