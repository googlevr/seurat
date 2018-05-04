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

#include "seurat/image/atlaser_util.h"

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/base/util.h"

namespace seurat {
namespace image {

using ion::math::Point2i;
using ion::math::Range2i;
using ion::math::Vector2i;

void ComputeAtlasLayout(absl::Span<const Vector2i> tile_sizes, int atlas_width,
                        int* height, absl::Span<Point2i> tile_origins) {
  DCHECK_EQ(tile_sizes.size(), tile_origins.size());
  const int num_tiles = tile_sizes.size();

  // Pair tile size with ID.
  std::vector<std::pair<int, Vector2i>> sorted_tile_sizes(num_tiles);
  for (int i = 0; i < num_tiles; ++i) {
    sorted_tile_sizes[i] = std::pair<int, Vector2i>(i, tile_sizes[i]);
  }
  // Sort the pairs, first by decreasing height then by decreasing width.
  std::sort(
      sorted_tile_sizes.begin(), sorted_tile_sizes.end(),
      [](const std::pair<int, Vector2i>& a, const std::pair<int, Vector2i>& b) {
        if (a.second[1] == b.second[1])
          return a.second[0] > b.second[0];
        else
          return a.second[1] > b.second[1];
      });

  Point2i next_free(0, 0);
  int row_height = 0;
  for (const auto& index_and_size : sorted_tile_sizes) {
    // Tiles must have non-zero width and height.
    CHECK_LT(0, index_and_size.second[0]) << "Tiles must have non-zero width.";
    CHECK_LT(0, index_and_size.second[1]) << "Tiles must have non-zero height.";
    if (next_free[0] + index_and_size.second[0] > atlas_width) {
      // The tile doesn't fit into the current row. Start a new row.
      next_free[0] = 0;
      next_free[1] += row_height;
    }
    if (next_free[0] == 0) {
      // The tile defines the height of the new row.
      row_height = index_and_size.second[1];
    }
    tile_origins[index_and_size.first] = next_free;
    next_free[0] += index_and_size.second[0];
  }
  *height = next_free[1] + row_height;
}

Vector2i ComputeAtlasSize(absl::Span<const Vector2i> tile_sizes,
                          absl::Span<const Point2i> tile_origins) {
  Range2i total_range;
  for (int i = 0; i < tile_sizes.size(); ++i) {
    total_range.ExtendByRange(
        Range2i::BuildWithSize(tile_origins[i], tile_sizes[i]));
  }
  return total_range.GetSize();
}

}  // namespace image
}  // namespace seurat
