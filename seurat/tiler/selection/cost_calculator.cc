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

#include "ion/base/logging.h"
#include "seurat/base/parallel.h"
#include "seurat/tiler/selection/selection_util.h"

namespace seurat {
namespace tiler {
namespace selection {

void CostCalculator::ComputeSolutionCost(const SelectionProblem &problem,
                                         absl::Span<const double> multipliers,
                                         absl::Span<const int> items,
                                         double *primal_cost, double *dual_cost,
                                         absl::Span<double> total_weight) {
  const int num_items = items.size();
  const int num_weights = problem.items->GetNumWeights();
  DCHECK_EQ(multipliers.size(), num_weights);
  DCHECK_EQ(total_weight.size(), num_weights);
  DCHECK_EQ(problem.capacity.size(), num_weights);

  // Allocate more chunks than threads to improve load balancing. Don't let the
  // number of chunks exceed the number of items, or the chunk processing will
  // go out of bounds.
  int num_chunks = std::min(thread_count_ * 2, num_items);
  // Round up.
  int chunk_size = (items.size() + num_chunks - 1) / num_chunks;
  // Neither num_chunks nor chunk_size can be chosen and fixed; they depend on
  // each other. Rounding up chunk_size sometimes results in the last chunk
  // overshooting the legal range of num_items, while rounding down sometimes
  // leads to missing part of the last block.
  //
  // Consider 5 items and 2 threads, so we desire four chunks. Rounding up the
  // chunk size leads to block size 2. Then the blocks are [0,1], [2,3], [4,5]
  // and [6,7] so the last block is out of bounds. Rounding down results in four
  // blocks of size 1, missing the last item.
  //
  // Therefore the chunk size rounding mode only influences the size of task for
  // each thread. Adjusting the number of chunks to match the block size
  // prevents overshoot and guarantees coverage.
  num_chunks = (num_items + chunk_size - 1) / chunk_size;

  partial_costs_.resize(num_chunks);

  // Process all chunks of items to compute their primal cost & total weight.
  base::BalancedParallelFor(thread_count_, num_chunks, [&](int chunk_id) {
    int start = chunk_size * chunk_id;
    int size = std::min(start + chunk_size, num_items) - start;
    absl::Span<const int> items_in_chunk = items.subspan(start, size);
    partial_costs_[chunk_id].primal_cost =
        ComputeTotalCost(*problem.items, items_in_chunk);
    partial_costs_[chunk_id].total_weight.resize(num_weights);
    ComputeTotalWeight(*problem.items, items_in_chunk,
                       absl::MakeSpan(partial_costs_[chunk_id].total_weight));
  });
  // Sum over all chunks.
  *primal_cost = 0.0;
  std::fill(total_weight.begin(), total_weight.end(), 0.0);
  for (const auto &partial : partial_costs_) {
    *primal_cost += partial.primal_cost;
    for (int w = 0; w < num_weights; ++w) {
      total_weight[w] += partial.total_weight[w];
    }
  }
  *dual_cost =
      ComputeDualCost(problem, multipliers, total_weight, *primal_cost);
}

}  // namespace selection
}  // namespace tiler
}  // namespace seurat
