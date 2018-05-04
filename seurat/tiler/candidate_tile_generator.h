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

#ifndef VR_SEURAT_TILER_CANDIDATE_TILE_GENERATOR_H_
#define VR_SEURAT_TILER_CANDIDATE_TILE_GENERATOR_H_

#include <functional>
#include <vector>

#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/tiler/build_partition.h"
#include "seurat/tiler/geometry_model.h"
#include "seurat/tiler/geometry_solver.h"
#include "seurat/tiler/partitioner_stage.h"
#include "seurat/tiler/point_set.h"
#include "seurat/tiler/subdivision.h"
#include "seurat/tiler/tile.h"
#include "seurat/tiler/tile_resolver.h"

namespace seurat {
namespace tiler {

// An interface for generating candidate tiles for points at various
// cells of a Subdivision.
//
// Implementations are not thread safe, unless otherwise specified.
class CandidateTileGenerator {
 public:
  // A set of candidate tiles fit to the points in a single subdivision cell.
  struct CandidateTiles {
    // The set of tiles, in no particular order.
    std::vector<Tile> tiles;

    // The geometric-reconstruction error of each tile.
    std::vector<float> costs;
  };

  virtual ~CandidateTileGenerator() = default;

  // Initializes the CandidateTileGenerator with the PointSet, enabling
  // precomputation of acceleration structures to be reused over multiple calls
  // to GenerateCandidatePartitionings.
  virtual void Init(const PointSet& point_set) = 0;

  // Given a PointSet and its organization in a pyramid, computes a set of
  // candidate tiles for each of the indicated cells of the pyramid.
  //
  // Note that |cell_indices| and |candidates_per_cell| are parallel arrays
  // with corresponding values.  That is, candidates_per_cell[i] will have
  // candidate tiles for the points in cell_indices[i].
  virtual void GenerateCandidateTiles(
      const PointSet& point_set, const Subdivision& subdivision,
      absl::Span<const int> cell_indices,
      absl::Span<std::vector<CandidateTiles>> candidates_per_cell) = 0;
};

// A CandidateTileGenerator which incrementally computes candidate tile
// clusterings up to a maximum size.
class ExhaustiveCandidateTileGenerator : public CandidateTileGenerator {
 public:
  ExhaustiveCandidateTileGenerator(
      int max_partitions, std::shared_ptr<PartitionerStage> child_partitioner,
      std::shared_ptr<GeometrySolver> geometry_solver,
      std::shared_ptr<TileResolver> tile_resolver)
      : max_partitions_(max_partitions),
        child_partitioner_(std::move(child_partitioner)),
        geometry_solver_(std::move(geometry_solver)),
        tile_resolver_(std::move(tile_resolver)) {}

  // CandidateTileGenerator implementation.
  void Init(const PointSet& point_set) override {
    geometry_solver_->Init(point_set);
    child_partitioner_->Init(point_set);
    tile_resolver_->Init(point_set);
  }
  void GenerateCandidateTiles(
      const PointSet& point_set, const Subdivision& subdivision,
      absl::Span<const int> cell_indices,
      absl::Span<std::vector<CandidateTiles>> candidates_per_cell) override;

 private:
  // Converts the given partitions into tiles.
  //
  // Returns the empty set upon failure.
  CandidateTiles ToCandidateTiles(
      const std::vector<BuildPartition>& partitions) const;

  // The maximum number of partitions to use in a single candidate partitioning.
  const int max_partitions_;

  // The partitioner used to cluster points.
  const std::shared_ptr<PartitionerStage> child_partitioner_;

  // Defines the model of the proxy-geometry for the partitions.
  const std::shared_ptr<GeometrySolver> geometry_solver_;

  // Converts BuildPartitions into Tiles.
  const std::shared_ptr<TileResolver> tile_resolver_;
};

// Executes a set of other CandidateTileGenerators in parallel.
class ParallelCandidateTileGenerator : public CandidateTileGenerator {
 public:
  // Constructs a CandidateTileGenerator.
  //
  // Note that the |generator_factory| must supply instances which can be
  // simultaneously used on separate threads, even if the instances themselves
  // are not thread-safe.
  ParallelCandidateTileGenerator(
      int thread_count,
      std::function<std::shared_ptr<CandidateTileGenerator>(void)>
          generator_factory);

  // CandidateTileGenerator implementation.
  void Init(const PointSet& point_set) override;
  void GenerateCandidateTiles(
      const PointSet& point_set, const Subdivision& subdivision,
      absl::Span<const int> cell_indices,
      absl::Span<std::vector<CandidateTiles>> candidates_per_cell) override;

 private:
  // The pool of CandidateTileGenerators to use on each thread.
  std::vector<std::shared_ptr<CandidateTileGenerator>> generators_;
};

}  // namespace tiler
}  // namespace seurat
#endif  // VR_SEURAT_TILER_CANDIDATE_TILE_GENERATOR_H_
