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

#include <vector>

#include "ion/math/vector.h"
#include "gtest/gtest.h"

namespace seurat {
namespace image {

using ion::math::Point2i;
using ion::math::Vector2i;

namespace {

TEST(AtlaserUtil, ComputeAtlasLayout) {
  const Vector2i kTileSize(10, 12);
  std::vector<Vector2i> tile_sizes = {kTileSize};
  std::vector<Point2i> tile_origins(tile_sizes.size());
  int atlas_height;
  ComputeAtlasLayout(tile_sizes, kTileSize[0], &atlas_height,
                     absl::MakeSpan(tile_origins));

  EXPECT_EQ(1, tile_origins.size());
  EXPECT_EQ(kTileSize[1], atlas_height);
  EXPECT_EQ(Point2i(0, 0), tile_origins[0]);
}

TEST(AtlaserUtil, ComputeAtlasSize) {
  const Vector2i kTileSize(3, 4);
  std::vector<Vector2i> tile_sizes = {kTileSize};
  std::vector<Point2i> tile_origins(tile_sizes.size());
  int atlas_height;
  ComputeAtlasLayout(tile_sizes, kTileSize[0], &atlas_height,
                     absl::MakeSpan(tile_origins));
  Vector2i atlas_size = ComputeAtlasSize(tile_sizes, tile_origins);
  EXPECT_EQ(kTileSize, atlas_size);
}

}  // namespace

}  // namespace image
}  // namespace seurat
