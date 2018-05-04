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

#include <vector>

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "gtest/gtest.h"
#include "seurat/base/color.h"
#include "seurat/base/util.h"
#include "seurat/image/fixed_width_atlaser.h"

using ion::math::Point2i;
using ion::math::Range2i;
using ion::math::Vector2i;
using seurat::base::Color3ui8;

namespace seurat {
namespace image {

namespace {

// Checks the sequence is strictly increasing. Note std::is_sorted is nearly
// sufficient but allows duplicates.
template <class ForwardIt>
bool is_increasing(ForwardIt first, ForwardIt last) {
  while (first != last) {
    ForwardIt next = first;
    ++next;
    if (next != last && !(*first < *next)) {
      return false;
    }
    first = next;
  }

  return true;
}

}  // namespace

// Test splitting up a tile set over multiple textures.
TEST(TextureAtlas, ConstrainedSortedTexSizeAtlas) {
  // These values are a subset of tile sizes output from the tiler while
  // processing a data set.
  const int kTestTextureTiles[119][2] = {
      {38, 38}, {38, 38},  {10, 22},  {25, 38}, {37, 38}, {37, 37}, {37, 38},
      {13, 16}, {37, 37},  {38, 38},  {38, 38}, {37, 38}, {15, 15}, {37, 38},
      {37, 37}, {11, 12},  {38, 38},  {38, 38}, {37, 37}, {38, 38}, {38, 38},
      {24, 37}, {38, 38},  {38, 38},  {38, 38}, {38, 38}, {38, 38}, {37, 38},
      {37, 38}, {37, 38},  {38, 38},  {38, 38}, {38, 38}, {37, 38}, {14, 16},
      {20, 23}, {20, 30},  {22, 38},  {17, 38}, {15, 16}, {26, 38}, {8, 9},
      {38, 38}, {20, 37},  {25, 39},  {5, 7},   {38, 38}, {37, 38}, {38, 38},
      {38, 38}, {37, 38},  {15, 15},  {37, 38}, {38, 38}, {26, 28}, {8, 15},
      {17, 24}, {31, 38},  {34, 38},  {37, 38}, {6, 8},   {7, 9},   {37, 38},
      {46, 48}, {43, 47},  {32, 42},  {43, 49}, {42, 55}, {20, 27}, {46, 47},
      {46, 58}, {42, 120}, {37, 53},  {51, 59}, {51, 55}, {35, 78}, {41, 119},
      {47, 51}, {48, 62},  {47, 85},  {50, 63}, {46, 49}, {37, 59}, {54, 56},
      {46, 49}, {44, 49},  {49, 51},  {41, 88}, {50, 54}, {51, 65}, {54, 92},
      {50, 69}, {48, 58},  {48, 53},  {49, 49}, {44, 49}, {46, 48}, {47, 48},
      {16, 34}, {35, 50},  {49, 57},  {46, 47}, {27, 50}, {46, 49}, {38, 85},
      {48, 51}, {19, 49},  {40, 54},  {50, 68}, {51, 64}, {50, 71}, {43, 56},
      {34, 48}, {46, 47},  {42, 114}, {18, 58}, {28, 71}, {47, 56}, {33, 75},
  };

  std::vector<Vector2i> tiles;
  for (const auto& size : kTestTextureTiles) {
    tiles.push_back(Vector2i(size[0], size[1]));
  }
  const Vector2i kAtlasSizeTarget(256, 128);
  std::vector<int> chart_splits =
      TextureAtlas::FindAtlasSplits(tiles, kAtlasSizeTarget);
  // No empty charts allowed, splits must be in order.
  EXPECT_EQ(is_increasing(chart_splits.begin(), chart_splits.end()), true);
  // All tiles must be included.
  EXPECT_EQ(chart_splits.back(), tiles.size());
  // The tile set can't fit on a single chart.
  EXPECT_GT(tiles.size(), 1);

  FixedWidthAtlaser atlaser(kAtlasSizeTarget);

  // All resulting charts must be under the texture size limit.
  int first_chart_tile = 0;
  for (int end_chart_tile : chart_splits) {
    absl::Span<const Vector2i> tile_subset = absl::MakeConstSpan(tiles).subspan(
        first_chart_tile, end_chart_tile - first_chart_tile);
    Vector2i subset_size;
    std::vector<Point2i> origins(tile_subset.size());
    atlaser.LayoutTiles(tile_subset, &subset_size, absl::MakeSpan(origins));
    EXPECT_LE(subset_size[0], kAtlasSizeTarget[0])
        << subset_size << " overflowed chart at tile " << end_chart_tile;
    EXPECT_LE(subset_size[1], kAtlasSizeTarget[1])
        << subset_size << " overflowed chart at tile " << end_chart_tile;

    first_chart_tile = end_chart_tile;
  }
}

}  // namespace image
}  // namespace seurat
