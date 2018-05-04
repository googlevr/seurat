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

#include "seurat/tiler/tile_weight_model.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/tiler/tile.h"

namespace seurat {
namespace tiler {
namespace {

using ion::math::Point3f;

// Makes a 0.2 x 0.2 quad at z = -1.
//
// The quad is offset to be centered at (0.5, 0.5, -1) such that
// perspective-foreshortening effects roughly even out.
Tile MakeSmallSquareTile() {
  const float kHalfSize = 0.1f;
  const Point3f offset(0.5f, 0.5f, 0.0f);
  Tile tile;
  tile.cell = 0;
  tile.quad = {{
      offset + Point3f(-kHalfSize, -kHalfSize, -1.0f),  //
      offset + Point3f(kHalfSize, -kHalfSize, -1.0f),   //
      offset + Point3f(kHalfSize, kHalfSize, -1.0f),    //
      offset + Point3f(-kHalfSize, kHalfSize, -1.0f)    //
  }};
  return tile;
}

// Like MakeSmallSquareTile(), but vertex order is reversed.
Tile MakeSmallSquareTileReversed() {
  const float kHalfSize = 0.1f;
  const Point3f offset(0.5f, 0.5f, 0.0f);
  Tile tile;
  tile.cell = 0;
  tile.quad = {{
      offset + Point3f(-kHalfSize, -kHalfSize, -1.0f),  //
      offset + Point3f(-kHalfSize, kHalfSize, -1.0f),   //
      offset + Point3f(kHalfSize, kHalfSize, -1.0f),    //
      offset + Point3f(kHalfSize, -kHalfSize, -1.0f)    //
  }};
  return tile;
}

// Returns a dense grid with z=-1 and with (x, y) in [-1.0, 1.0]^2.
//
// This is essentially one face of a cubemap.
std::vector<Tile> MakeTessellatedQuad() {
  const int kResolution = 10;
  std::vector<Tile> tiles;
  for (int y = 0; y < kResolution; ++y) {
    for (int x = 0; x < kResolution; ++x) {
      std::array<Point3f, 4> quad = {{Point3f(x, y, 0.0f),          //
                                      Point3f(x + 1, y, 0.0f),      //
                                      Point3f(x + 1, y + 1, 0.0f),  //
                                      Point3f(x, y + 1, 0.0f)}};
      for (Point3f& p : quad) {
        p /= static_cast<float>(kResolution);
        p -= {0.5f, 0.5f, 0.0f};
        p *= 2.0f;
        p[2] = -1.0f;
      }
      Tile tile;
      tile.cell = 0;
      tile.quad = quad;
      tiles.push_back(tile);
    }
  }

  return tiles;
}

TEST(TileWeightModelTest, TriangleWeightModel) {
  TriangleCountTileWeightModel model;

  std::vector<float> weight(1);
  EXPECT_TRUE(
      model.GetTileWeight(MakeSmallSquareTile(), absl::MakeSpan(weight)));

  // 2 triangles per tile.
  EXPECT_EQ(2.0f, weight[0]);
}

TEST(TileWeightModelTest, ProjectedAreaWeightModel) {
  ProjectedAreaTileWeightModel model;

  std::vector<float> weight(1);
  EXPECT_TRUE(
      model.GetTileWeight(MakeSmallSquareTile(), absl::MakeSpan(weight)));

  // This is a 0.1 by 0.1 square on a single face of a cubemap.
  const float kExpectedSize = 0.01 / 6.0f;
  EXPECT_NEAR(kExpectedSize, weight[0], 1e-4f);

  model.GetTileWeight(MakeSmallSquareTileReversed(), absl::MakeSpan(weight));
  EXPECT_NEAR(kExpectedSize, weight[0], 1e-4f);
}

TEST(TileWeightModelTest, CombinedWeightModel) {
  std::vector<std::unique_ptr<TileWeightModel>> models;
  models.emplace_back(new TriangleCountTileWeightModel);
  models.emplace_back(new ProjectedAreaTileWeightModel);
  CombinedTileWeightModel combined_model(std::move(models));

  EXPECT_EQ(2, combined_model.GetDimension());

  std::array<float, 2> weight;
  EXPECT_TRUE(combined_model.GetTileWeight(MakeSmallSquareTile(),
                                           absl::MakeSpan(weight)));

  EXPECT_EQ(2.0f, weight[0]);
  const float kExpectedSize = 0.01 / 6.0f;
  EXPECT_NEAR(kExpectedSize, weight[1], 1e-3f);
}

std::vector<Tile> MakeCubemapTiles() {
  std::vector<Tile> tiles;
  // Process tiles which form all 6 faces of a cube map.
  //
  // The total weight from all samples over all tiles should be near 1.0 for all
  // directions.
  for (int i = 0; i < 6; ++i) {
    int major_axis = i % 3;
    int sign = (i / 3) * 2 - 1;

    // A tile which is one face of an origin-centered cube with
    // half-side-length=1.
    Tile tile;
    std::fill(tile.quad.begin(), tile.quad.end(), Point3f::Zero());
    tile.quad[0][major_axis] = sign;
    tile.quad[0][(major_axis + 1) % 3] = -sign;
    tile.quad[0][(major_axis + 2) % 3] = -sign;

    tile.quad[1][major_axis] = sign;
    tile.quad[1][(major_axis + 1) % 3] = sign;
    tile.quad[1][(major_axis + 2) % 3] = -sign;

    tile.quad[2][major_axis] = sign;
    tile.quad[2][(major_axis + 1) % 3] = sign;
    tile.quad[2][(major_axis + 2) % 3] = sign;

    tile.quad[3][major_axis] = sign;
    tile.quad[3][(major_axis + 1) % 3] = -sign;
    tile.quad[3][(major_axis + 2) % 3] = sign;

    tiles.push_back(tile);
  }
  return tiles;
}

TEST(TileWeightModelTest, DirectionalOverdrawTileWeightModel_Cubemap) {
  const int kSamples = 20;
  const int kFieldOfViewRadians = 90;
  // The geometry is positioned around 1 unit from the origin, so scale the
  // headbox appropriately.
  const float kHeadboxRadius = 0.01f;
  std::unique_ptr<TileWeightModel> model =
      DirectionalOverdrawTileWeightModel::Build(kSamples, kFieldOfViewRadians,
                                                kHeadboxRadius);

  std::vector<float> total_weight(model->GetDimension(), 0.0f);
  std::vector<float> cur_weight(model->GetDimension());

  std::vector<Tile> tiles = MakeCubemapTiles();

  // The total weight from all samples over all tiles should be near 1.0 for all
  // directions.
  for (const auto& tile : tiles) {
    EXPECT_TRUE(model->GetTileWeight(tile, absl::MakeSpan(cur_weight)));
    for (int w = 0; w < model->GetDimension(); ++w) {
      total_weight[w] += cur_weight[w];
    }
  }

  for (int w = 0; w < model->GetDimension(); ++w) {
    // Total overdraw for all tiles should be near 1 for each view.
    //
    // Since this is very-much a heuristic, allow any values between 0.95 and
    // 1.05.
    EXPECT_NEAR(1.0f, total_weight[w], 0.05f);
  }
}

TEST(TileWeightModelTest, DirectionalOverdrawTileWeightModel_SmallTile) {
  const int kSamples = 20;
  const int kFieldOfViewRadians = 90;
  const float kHeadboxRadius = 0.0f;
  std::unique_ptr<TileWeightModel> model =
      DirectionalOverdrawTileWeightModel::Build(kSamples, kFieldOfViewRadians,
                                                kHeadboxRadius);

  std::vector<float> weight(model->GetDimension());

  // Make a *really* small tile.
  const float kHalfSize = 0.001f;
  const Point3f offset(0.5f, 0.5f, 0.0f);
  Tile tile;
  tile.cell = 0;
  tile.quad = {{
      offset + Point3f(-kHalfSize, -kHalfSize, -1.0f),  //
      offset + Point3f(kHalfSize, -kHalfSize, -1.0f),   //
      offset + Point3f(kHalfSize, kHalfSize, -1.0f),    //
      offset + Point3f(-kHalfSize, kHalfSize, -1.0f)    //
  }};

  EXPECT_TRUE(model->GetTileWeight(tile, absl::MakeSpan(weight)));

  float max_weight = *std::max_element(weight.begin(), weight.end());
  // At least one element should be non-zero.
  EXPECT_LT(0.0f, max_weight);
}

}  // namespace
}  // namespace tiler
}  // namespace seurat
