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

#include "seurat/tiler/subdivision_util.h"

#include "absl/types/span.h"

namespace seurat {
namespace tiler {

void GetCellsInDepthRange(const Subdivision& tree, int min_depth, int max_depth,
                          std::vector<int>* relevant_cells) {
  relevant_cells->clear();
  struct WorkItem {
    int node;
    int depth;
  };
  std::vector<WorkItem> to_visit;
  std::vector<int> child_nodes;

  tree.GetRoots(&child_nodes);
  to_visit.reserve(child_nodes.size());
  for (const int& root : child_nodes) {
    to_visit.push_back({root, 0});
  }

  while (!to_visit.empty()) {
    child_nodes.clear();

    const WorkItem item = to_visit.back();
    to_visit.pop_back();

    if (item.depth >= min_depth) {
      relevant_cells->push_back(item.node);
    }

    if (item.depth < max_depth) {
      tree.GetChildren(item.node, &child_nodes);
      for (int child : child_nodes) {
        to_visit.push_back({child, item.depth + 1});
      }
    }
  }
}

}  // namespace tiler
}  // namespace seurat
