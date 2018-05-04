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

#ifndef VR_SEURAT_TILER_SELECTION_TILER_H_
#define VR_SEURAT_TILER_SELECTION_TILER_H_

#include <memory>
#include <vector>

#include "seurat/tiler/candidate_tile_generator.h"
#include "seurat/tiler/point_set.h"
#include "seurat/tiler/selection/selection_solver.h"
#include "seurat/tiler/subdivision.h"
#include "seurat/tiler/tile_weight_model.h"
#include "seurat/tiler/tiler.h"

namespace seurat {
namespace tiler {

// Considers all cells of a hierarchical Subdivision (beyond some minimum depth)
// to:
//  1. Generate a set of candidate tiles (using a CandidateTileGenerator).
//  2. Select the best candidates to use, by encoding the selection-problem as a
//     SelectionProblem.
class SelectionTiler : public Tiler {
 public:
  SelectionTiler(int min_subdivision_level,
                 std::shared_ptr<Subdivision> subdivision,
                 std::shared_ptr<CandidateTileGenerator> candidate_generator,
                 std::shared_ptr<selection::SelectionSolver> selection_solver,
                 std::shared_ptr<TileWeightModel> tile_weight_model,
                 std::vector<float> max_weight, int thread_count)
      : min_subdivision_level_(min_subdivision_level),
        subdivision_(std::move(subdivision)),
        candidate_generator_(std::move(candidate_generator)),
        selection_solver_(std::move(selection_solver)),
        tile_weight_model_(std::move(tile_weight_model)),
        max_weight_(std::move(max_weight)),
        thread_count_(thread_count) {}
  ~SelectionTiler() override = default;

  std::vector<Tile> Run(const PointSet& points) override;

 private:
  // The minimum subdivision level.
  const int min_subdivision_level_;

  // The Subdivision to use for tiling.
  const std::shared_ptr<Subdivision> subdivision_;

  // Generates candidate tiles for each cell of the |subdivision_|.
  const std::shared_ptr<CandidateTileGenerator> candidate_generator_;

  // Used to select the best geometry to use.
  const std::shared_ptr<selection::SelectionSolver> selection_solver_;

  // Determines the weight of all tiles.
  const std::shared_ptr<TileWeightModel> tile_weight_model_;

  // The maximum total weight of all output tiles.
  const std::vector<float> max_weight_;

  // The number of threads to use.
  const int thread_count_;
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SELECTION_TILER_H_
