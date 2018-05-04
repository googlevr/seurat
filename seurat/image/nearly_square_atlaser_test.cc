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

#include <random>
#include <vector>

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "gtest/gtest.h"
#include "seurat/base/color.h"
#include "seurat/base/util.h"

namespace seurat {
namespace image {

using ion::math::Point2i;
using ion::math::Range2i;
using ion::math::Vector2i;

namespace {

TEST(NearlySquareAtlaser, GetAtlasSizeTarget) {
  NearlySquareAtlaser atlaser;
  EXPECT_EQ(Vector2i::Zero(), atlaser.GetAtlasSizeTarget());
}

TEST(NearlySquareAtlaser, LayoutSingleTile) {
  std::vector<Vector2i> tile_sizes;
  tile_sizes.push_back({10, 10});

  Vector2i size;
  std::vector<Point2i> tile_origins(tile_sizes.size());
  NearlySquareAtlaser atlaser;
  atlaser.LayoutTiles(tile_sizes, &size, absl::MakeSpan(tile_origins));

  EXPECT_EQ(1, tile_origins.size());
  EXPECT_EQ(Vector2i(10, 10), size);
  EXPECT_EQ(Point2i(0, 0), tile_origins[0]);
}

TEST(NearlySquareAtlaser, LayoutTiles) {
  const int kNumRuns = 100;
  std::default_random_engine engine(12345);
  std::uniform_int_distribution<int> dist_size(1, 32);
  std::uniform_int_distribution<int> dist_num_tiles(1, 64);

  for (int run = 0; run < kNumRuns; ++run) {
    // Create a set of random tiles (we just need the sizes).
    std::vector<Vector2i> tile_sizes;
    for (int i = 0; i < dist_num_tiles(engine); ++i) {
      tile_sizes.push_back({dist_size(engine), dist_size(engine)});
    }

    // Layout the tiles.
    Vector2i size;
    std::vector<Point2i> tile_origins(tile_sizes.size());
    NearlySquareAtlaser atlaser;
    atlaser.LayoutTiles(tile_sizes, &size, absl::MakeSpan(tile_origins));
    EXPECT_EQ(tile_sizes.size(), tile_origins.size());

    // Check that tiles do not overlap.
    for (int i = 0; i < tile_sizes.size(); ++i) {
      Range2i tile_i_range = Range2i::BuildWithSize(
          tile_origins[i], tile_sizes[i] - Vector2i(1, 1));
      for (int j = i + 1; j < tile_sizes.size(); ++j) {
        Range2i tile_j_range = Range2i::BuildWithSize(
            tile_origins[j], tile_sizes[j] - Vector2i(1, 1));
        EXPECT_FALSE(tile_i_range.IntersectsRange(tile_j_range));
      }
    }

    // Verify that the size of the resulting atlas is equal to the specified
    // size.
    Range2i total_size;
    for (int i = 0; i < tile_sizes.size(); ++i) {
      total_size.ExtendByRange(
          Range2i::BuildWithSize(tile_origins[i], tile_sizes[i]));
    }
    EXPECT_EQ(Point2i::Zero(), total_size.GetMinPoint());
    EXPECT_EQ(total_size.GetSize(), size);
  }
}

}  // namespace

}  // namespace image
}  // namespace seurat
