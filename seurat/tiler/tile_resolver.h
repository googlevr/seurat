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

#ifndef VR_SEURAT_TILER_TILE_RESOLVER_H_
#define VR_SEURAT_TILER_TILE_RESOLVER_H_

#include "seurat/tiler/build_partition.h"
#include "seurat/tiler/geometry_model.h"
#include "seurat/tiler/point_set.h"
#include "seurat/tiler/subdivision.h"
#include "seurat/tiler/tile.h"

namespace seurat {
namespace tiler {

// An interface for converting a BuildPartition into a Tile.
//
// This converts from the *implicit-geometry* of a partition into
// *explicit-geometry* of a tile.
class TileResolver {
 public:
  virtual ~TileResolver() = default;

  // This must always be called before Resolve().  TileResolvers may use this
  // to precompute acceleration-structures based on the points.
  virtual void Init(const PointSet& point_set) {
    // Do nothing by default.
  }

  virtual bool Resolve(const BuildPartition& partition, Tile* tile) const = 0;
};

// Generates tile quads by intersecting the partition's plane with its
// subdivision cell's "rails".
class RailTileResolver : public TileResolver {
 public:
  explicit RailTileResolver(std::shared_ptr<Subdivision> subdivision)
      : subdivision_(std::move(subdivision)) {}

  // TileResolver implementation.
  ~RailTileResolver() override = default;
  void Init(const PointSet& point_set) override;
  bool Resolve(const BuildPartition& partition, Tile* tile) const override;

 private:
  const std::shared_ptr<Subdivision> subdivision_;
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_TILE_RESOLVER_H_
