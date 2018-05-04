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

#include "seurat/tiler/subset_geometry_solver.h"

#include <algorithm>
#include <numeric>
#include <vector>

#include "seurat/base/util.h"

namespace seurat {
namespace tiler {

bool SubsetGeometrySolver::FitModel(absl::Span<const int> point_indices,
                                    GeometryModel *model) const {
  const int num_points = point_indices.size();
  std::vector<int> subset(std::min(max_point_count_, num_points));
  float scale = num_points / static_cast<float>(subset.size());
  for (int i = 0; i < subset.size(); ++i) {
    subset[i] = point_indices[i * scale];
  }
  return delegate_->FitModel(subset, model);
}

}  // namespace tiler
}  // namespace seurat
