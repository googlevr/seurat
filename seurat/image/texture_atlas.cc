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

#include "seurat/image/texture_atlas.h"

#include <algorithm>
#include <numeric>
#include <vector>

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/base/util.h"
#include "seurat/image/fixed_width_atlaser.h"

namespace seurat {
namespace image {

using ion::math::Point2i;
using ion::math::Range2i;
using ion::math::Vector2i;

std::vector<int> TextureAtlas::FindAtlasSplits(
    absl::Span<const Vector2i> tile_sizes, const Vector2i& max_size) {
  int start_tile = 0;
  std::vector<int> atlas_splits;
  CHECK_GT(max_size[0], 0);
  CHECK_GT(max_size[1], 0);

  // Scratch space to store each atlasing result.
  std::vector<Point2i> tile_origins(tile_sizes.size());

  // Build a help vector of indices to facilitate std::upper_bound.
  std::vector<int> tile_indices(tile_sizes.size());
  std::iota(tile_indices.begin(), tile_indices.end(), 0);

  // Search process:
  //
  // 1. Initialize start_tile to zero.
  // 2. Loop until start_tile reaches the count of tiles, packing as many tiles
  //    as possible into each chart.
  //    a. Start with the bracket [start_tile, start_tile+1).
  //    b. Quadratic probe up (double tile count) until the chart overflows.
  //    -  Now, we have a bracket of the maximal tile, somewhere between the
  //       tile count prior to doubling and the current tile count that
  //       overflowed the chart.
  //    c. Binary search the bracket until we find the last tile that we can add
  //       without overflowing the chart.
  //    d. Output the interval [start_tile, last_tile).
  //    e. Advance start_tile to last_tile.
  //
  // [ start tile
  // | first tile
  // * maximal tile
  //
  //            start                        maximal
  // 0----------[)---------------------------*-------------------N
  // 0----------[-)--------------------------*-------------------N
  // 0----------[---)------------------------*-------------------N
  // 0----------[-------)--------------------*-------------------N
  // 0----------[---------------)------------*-------------------N
  // 0----------[---------------|------------*--)----------------N
  // 0----------[-----------------------|----*--)----------------N
  // 0----------[---------------------------|*--)----------------N
  // 0----------[----------------------------*|-)----------------N
  //
  do {
    // First tile in the bracket.
    int first_tile = start_tile;
    // Last tile in the bracket.
    int last_tile = start_tile + 1;
    // Prior state of the last tile, used to update the bracket start after the
    // quadratic probing.
    int prior_last = last_tile;

    // Find the end of the bracket with quadratic probing.
    Vector2i subset_size;
    FixedWidthAtlaser atlaser(max_size);
    while (last_tile < tile_sizes.size() && subset_size[0] <= max_size[0] &&
           subset_size[1] <= max_size[1]) {
      // Quadratic probe up past current last tile.
      prior_last = last_tile;
      last_tile = first_tile + (last_tile - first_tile) * 2;
      last_tile = std::min(static_cast<int>(tile_sizes.size()), last_tile);

      absl::Span<const Vector2i> tile_subset =
          absl::MakeConstSpan(tile_sizes)
              .subspan(start_tile, last_tile - start_tile);
      tile_origins.resize(tile_subset.size());
      atlaser.LayoutTiles(tile_subset, &subset_size,
                          absl::MakeSpan(tile_origins));
    }

    // Adjust the start of the bracket. The last range didn't overflow the
    // chart, so advance the start of the search to it.
    first_tile = prior_last;

    // The binary search compares a tile index in the search range with the
    // maximum chart size, by taking the size of the range of tiles from the
    // start tile to the tile under test. This implicitly partitions the tiles
    // into two sets: one where the tiles fit in the chart - the overflow lambda
    // returns false - and one where the tiles overflow the chart as indicated
    // by the lambda returning true. Find the first tile that overflows.
    auto maximal_tile_location = std::upper_bound(
        tile_indices.begin() + first_tile, tile_indices.begin() + last_tile,
        max_size,
        [tile_sizes, start_tile, &subset_size, &tile_origins](
            const Vector2i& max_tileset_size, int test_tile_index) {
          // Note this subset is inclusive of |test_tile_index|, which is subtly
          // different from the exclusive treatment of the resulting maximal
          // tile index, as the maximal tile overflows the chart. This internal
          // test checks if a particular tile included in the range overflows
          // the chart.
          absl::Span<const Vector2i> tile_subset =
              absl::MakeConstSpan(tile_sizes)
                  .subspan(start_tile, test_tile_index - start_tile + 1);
          tile_origins.resize(tile_subset.size());
          FixedWidthAtlaser atlaser(max_tileset_size);
          atlaser.LayoutTiles(tile_subset, &subset_size,
                              absl::MakeSpan(tile_origins));
          bool chart_too_small = max_tileset_size[0] < subset_size[0] ||
                                 max_tileset_size[1] < subset_size[1];
          return chart_too_small;
        });
    int maximal_tile = (maximal_tile_location == tile_indices.end())
                           ? tile_indices.size()
                           : *maximal_tile_location;

    // Output the split point.
    atlas_splits.push_back(maximal_tile);

    // Advance past the consumed tiles.
    start_tile = maximal_tile;
  } while (start_tile < tile_sizes.size());

  return atlas_splits;
}

}  // namespace image
}  // namespace seurat
