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

#include "seurat/tiler/candidate_tile_generator.h"

#include <algorithm>
#include <atomic>

#include "seurat/base/parallel.h"
#include "seurat/base/progress.h"
#include "seurat/tiler/partitioner_stages.h"

namespace seurat {
namespace tiler {

using ion::math::Point3f;

CandidateTileGenerator::CandidateTiles
ExhaustiveCandidateTileGenerator::ToCandidateTiles(
    const std::vector<BuildPartition>& partitions) const {
  std::vector<Tile> tiles;
  std::vector<float> cost;
  tiles.reserve(partitions.size());

  for (const auto& bp : partitions) {
    if (bp.GetPointIndices().empty()) {
      // Ignore empty partitions.
      continue;
    }
    Tile tile;
    if (!tile_resolver_->Resolve(bp, &tile)) {
      return {};
    }
    tiles.push_back(tile);
    cost.push_back(bp.GetTotalError());
  }

  return CandidateTiles{std::move(tiles), std::move(cost)};
}

void ExhaustiveCandidateTileGenerator::GenerateCandidateTiles(
    const PointSet& point_set, const Subdivision& subdivision,
    absl::Span<const int> cell_indices,
    absl::Span<std::vector<CandidateTiles>> candidates_per_cell) {
  tile_resolver_->Init(point_set);
  child_partitioner_->Init(point_set);
  for (int i = 0; i < cell_indices.size(); ++i) {
    int cell = cell_indices[i];
    std::vector<CandidateTiles>& candidate_tiles = candidates_per_cell[i];
    candidate_tiles.clear();

    absl::Span<const int> points_to_partition =
        subdivision.GetPointsInCell(cell);
    if (points_to_partition.empty()) {
      // This cell has no points.  Add the empty set as a valid set of candidate
      // tiles.
      CandidateTiles empty_candidate;
      candidate_tiles.push_back(empty_candidate);
      continue;
    }

    std::vector<BuildPartition> partitions;
    GeometryModel base_model;
    base_model.cell = cell;
    partitions.push_back(BuildPartition(base_model));
    for (int point : points_to_partition) {
      partitions.back().AddPoint(point, std::numeric_limits<float>::infinity());
    }

    for (int partition_count = 1; partition_count <= max_partitions_;
         ++partition_count) {
      child_partitioner_->Run(point_set, absl::MakeSpan(partitions));

      // The child partitioner must not add or remove partitions.
      //
      // If this fails, it means the ExhaustiveCandidateTileGenerator was
      // *misconfigured* with a child_partitioner_factory returning a
      // PartitionerStage which does not satisfy this requirement.
      CHECK_EQ(partition_count, partitions.size())
          << "ExhaustiveCandidateTileGenerator must be configured to wrap "
             "a child PartitionerStage which uses all partitions passed to it, "
             "but never adds or removes partitions";

      CandidateTiles candidates = ToCandidateTiles(partitions);
      if (!candidates.tiles.empty()) {
        candidate_tiles.push_back(std::move(candidates));
      }

      partitions.push_back(BuildPartition(base_model));
    }
  }
}

ParallelCandidateTileGenerator::ParallelCandidateTileGenerator(
    int thread_count,
    std::function<std::shared_ptr<CandidateTileGenerator>(void)>
        generator_factory)
    : generators_(thread_count) {
  std::generate(generators_.begin(), generators_.end(), generator_factory);
}

void ParallelCandidateTileGenerator::Init(const PointSet& point_set) {
  for (auto& generator : generators_) {
    generator->Init(point_set);
  }
}

void ParallelCandidateTileGenerator::GenerateCandidateTiles(
    const PointSet& point_set, const Subdivision& subdivision,
    absl::Span<const int> cell_indices,
    absl::Span<std::vector<CandidateTiles>> candidates_per_cell) {
  // Parallelize by treating cell_indices as a stack, with head at
  // next_cell_index.
  //
  // Each thread pops off the stack & processes that cell into its corresponding
  // output vector in candidate_partitionings.
  base::ScopedProgressRange progress("Generating mesh", cell_indices.size());
  std::atomic<int> next_cell_index(0);
  base::ParallelFor(generators_.size(), generators_.size(), [&](int tid) {
    CandidateTileGenerator& generator = *generators_[tid];
    while (true) {
      int current_cell_index = next_cell_index.fetch_add(1);

      if (current_cell_index >= cell_indices.size()) {
        return;
      }

      generator.GenerateCandidateTiles(
          point_set, subdivision, cell_indices.subspan(current_cell_index, 1),
          candidates_per_cell.subspan(current_cell_index, 1));
      progress.IncrementRange(1);
    }
  });
}

}  // namespace tiler
}  // namespace seurat
