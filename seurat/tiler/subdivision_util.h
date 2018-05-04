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

#ifndef VR_SEURAT_TILER_SUBDIVISION_UTIL_H_
#define VR_SEURAT_TILER_SUBDIVISION_UTIL_H_

#include <utility>
#include <vector>

#include "absl/types/span.h"
#include "seurat/tiler/point_set.h"
#include "seurat/tiler/subdivision.h"

namespace seurat {
namespace tiler {

// Returns all cells of the |tree| with
//   |min_depth| <= depth <= |max_depth|
//
// Returned cells are ordered according to a topological sort (top to bottom).
void GetCellsInDepthRange(const Subdivision& tree, int min_depth, int max_depth,
                          std::vector<int>* relevant_cells);

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SUBDIVISION_UTIL_H_
