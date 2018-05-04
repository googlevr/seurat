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

#ifndef VR_SEURAT_TILER_TILE_WEIGHT_MODEL_H_
#define VR_SEURAT_TILER_TILE_WEIGHT_MODEL_H_

#include <vector>

#include "absl/types/span.h"
#include "seurat/tiler/build_partition.h"
#include "seurat/tiler/geometry_model.h"
#include "seurat/tiler/point_set.h"
#include "seurat/tiler/tile.h"

namespace seurat {
namespace tiler {

// Models the "weight" of a tile.  This is a multidimensional quantity
// representing the cost of rendering, compressing, etc.
//
// The term "weight" is used here (as opposed to "cost") to distinguish
// from the geometric reconstruction cost typically associated with a set of
// tiles.
class TileWeightModel {
 public:
  virtual ~TileWeightModel() = default;

  // The number of components to the weight.
  virtual int GetDimension() const = 0;

  // Computes and returns the |weight| of the given |tile|.
  //
  // The given |weight| array must be of the correct dimension.
  //
  // Returns true upon success, false if the tile is ill-formed and should be
  // discarded.
  virtual bool GetTileWeight(const Tile& tile,
                             absl::Span<float> weight) const = 0;
};

// Weighs tiles based on the triangle-count required to render them.
class TriangleCountTileWeightModel : public TileWeightModel {
 public:
  ~TriangleCountTileWeightModel() override = default;

  int GetDimension() const override { return 1; }
  bool GetTileWeight(const Tile& tile,
                     absl::Span<float> weight) const override {
    DCHECK_EQ(weight.size(), GetDimension());
    // Two triangles per tile.
    weight[0] = 2.0f;
    return true;
  }
};

// Weighs tiles based on their area when projected onto the unit sphere.
//
// A weight of 1.0 indicates a tile which covers the unit-sphere.
class ProjectedAreaTileWeightModel : public TileWeightModel {
 public:
  ~ProjectedAreaTileWeightModel() override = default;

  int GetDimension() const override { return 1; }
  bool GetTileWeight(const Tile& tile, absl::Span<float> weight) const override;
};

// Concatenates multiple models.
class CombinedTileWeightModel : public TileWeightModel {
 public:
  CombinedTileWeightModel(std::vector<std::unique_ptr<TileWeightModel>> models)
      : models_(std::move(models)) {}
  ~CombinedTileWeightModel() override = default;

  int GetDimension() const override;
  bool GetTileWeight(const Tile& tile, absl::Span<float> weight) const override;

 private:
  std::vector<std::unique_ptr<TileWeightModel>> models_;
};

// Models the overdraw/fillrate required to render a set of tiles from different
// poses.
class DirectionalOverdrawTileWeightModel : public TileWeightModel {
 public:
  // Builds a Model which uses the specified number of directional |samples|,
  // with each sample modeling the overdraw as measured from a camera with the
  // specified field_of_view.
  static std::unique_ptr<TileWeightModel> Build(int samples,
                                                float field_of_view,
                                                float headbox_radius);

  ~DirectionalOverdrawTileWeightModel() override = default;

  int GetDimension() const override { return directions_.size(); }
  bool GetTileWeight(const Tile& tile, absl::Span<float> weight) const override;

 private:
  DirectionalOverdrawTileWeightModel(
      std::vector<ion::math::Vector3f> directions, float field_of_view_radians,
      float headbox_radius)
      : directions_(std::move(directions)),
        cos_half_field_of_view_(std::cos(field_of_view_radians / 2.0f)),
        headbox_radius_(headbox_radius) {}

  const std::vector<ion::math::Vector3f> directions_;
  const float cos_half_field_of_view_;
  // The radius of the (spherical) viewing volume.
  const float headbox_radius_;
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_TILE_WEIGHT_MODEL_H_
