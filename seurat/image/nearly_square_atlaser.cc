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

#include "seurat/image/nearly_square_atlaser.h"

#include <algorithm>
#include <numeric>

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/base/util.h"
#include "seurat/image/atlaser_util.h"

namespace seurat {
namespace image {

using ion::math::Point2i;
using ion::math::Range2i;
using ion::math::Vector2i;

namespace {

// Computes the required width for an atlas that can hold the tiles with the
// given |tile_sizes| such that the result will be approximately-square.
int CalculateAtlasWidth(absl::Span<const Vector2i> tile_sizes) {
  if (tile_sizes.empty()) return 0;
  int total_area = 0;
  for (const auto& tile_size : tile_sizes) {
    total_area += tile_size[0] * tile_size[1];
  }
  // Compute a width that will result in an approximately square atlas.  The
  // actual length will be different, because of overhead. Make the atlas at
  // least as wide as the widest tile.
  const Vector2i widest_tile = *std::max_element(
      tile_sizes.begin(), tile_sizes.end(),
      [](const Vector2i& a, const Vector2i& b) { return a[0] < b[0]; });
  const int sqrt_area = static_cast<int>(std::ceil(std::sqrt(total_area)));
  const int width = std::max(widest_tile[0], sqrt_area);
  return width;
}

}  // namespace

void NearlySquareAtlaser::LayoutTiles(
    absl::Span<const ion::math::Vector2i> tile_sizes,
    ion::math::Vector2i* total_size,
    absl::Span<ion::math::Point2i> tile_origins) const {
  DCHECK_EQ(tile_sizes.size(), tile_origins.size());
  const int initial_width_estimate = CalculateAtlasWidth(tile_sizes);
  int height;
  ComputeAtlasLayout(tile_sizes, initial_width_estimate, &height, tile_origins);
  *total_size = ComputeAtlasSize(tile_sizes, tile_origins);
}

}  // namespace image
}  // namespace seurat
