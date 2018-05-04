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

#ifndef VR_SEURAT_TILER_SELECTION_COST_CALCULATOR_H_
#define VR_SEURAT_TILER_SELECTION_COST_CALCULATOR_H_

#include "absl/types/span.h"
#include "seurat/tiler/selection/selection_problem.h"

namespace seurat {
namespace tiler {
namespace selection {

// Computes aggregate values for solutions to SelectionProblems.
//
// For efficiency, allocations for scratch space are cached.
//
// This is not thread safe.
class CostCalculator {
 public:
  explicit CostCalculator(int thread_count) : thread_count_(thread_count) {}

  // Given a problem and solution |multipliers| and selected |items|, computes
  // the primal & dual costs and the total weight of the solution.
  void ComputeSolutionCost(const SelectionProblem& problem,
                           absl::Span<const double> multipliers,
                           absl::Span<const int> items, double* primal_cost,
                           double* dual_cost, absl::Span<double> total_weight);

 private:
  // Partially computed values for some subset of all items in a solution.
  struct PartialCost {
    std::vector<double> total_weight;
    double primal_cost;
  };

  // The maximum number of threads to use.
  const int thread_count_;

  // Total costs computed for subsets of the current solution being processed.
  //
  // These are computed in parallel.
  std::vector<PartialCost> partial_costs_;
};

}  // namespace selection
}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SELECTION_COST_CALCULATOR_H_
