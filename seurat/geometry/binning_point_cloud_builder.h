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

#ifndef VR_SEURAT_GEOMETRY_BINNING_POINT_CLOUD_BUILDER_H_
#define VR_SEURAT_GEOMETRY_BINNING_POINT_CLOUD_BUILDER_H_

#include <array>
#include <deque>
#include <unordered_map>
#include <vector>

#include "ion/math/matrix.h"
#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "seurat/base/projective_camera.h"
#include "seurat/geometry/point_cloud_builder.h"

namespace seurat {
namespace geometry {

// This point cloud builder merges point sets by binning the points into six
// perspective grids, one for each cube-face camera at the origin. If a new
// point maps to a grid cell that is already occupied, that point will be
// discarded.
class BinningPointCloudBuilder : public PointCloudBuilder {
 public:
  BinningPointCloudBuilder(int thread_count,
                           const ion::math::Vector3i& resolution,
                           float near_clip);
  ~BinningPointCloudBuilder() override = default;

  // Computes the resolution of the binning grid and its near clipping plane for
  // the given |headbox| and |pixels_per_degree| target resolution.
  static void ComputeResolutionAndNearClip(const ion::math::Range3f& headbox,
                                           float pixels_per_degree,
                                           ion::math::Vector3i* resolution,
                                           float* near_clip);

  // PointCloudBuilder implementation.
  base::Status AddPoints(
      const std::vector<ion::math::Point3f>& positions) override;
  base::Status AddPointsWithWeights(
      const std::vector<ion::math::Point3f>& positions,
      const std::vector<float>& weights) override;
  void GetPositionsAndWeights(std::vector<ion::math::Point3f>* positions,
                              std::vector<float>* weights) override;

 private:
  // Data stored in a single cell (bin) of the six perspective grids.
  struct CellEntry {
    CellEntry(ion::math::Point3f position, float weight)
        : position(position), weight(weight) {}

    // Position of the representative point.
    ion::math::Point3f position;

    // Weight of the point. This is the sum of the weights of the points that
    // were binned into this cell.
    float weight;
  };

  // Computes the 64-bit index of the cell that contains |position_world|. The
  // index is a globally unique value across all six perspective grids. An error
  // status is returned, if an error occurred. Otherwise, the cell index is
  // stored in |cell_index|.
  base::Status CellIndexFromPoint(const ion::math::Point3f& position_world,
                                  int64* cell_index) const;

  // Number of threads used.
  const int thread_count_;

  // Resolution (x, y, z) of each of the six perspective grids.
  const ion::math::Vector3i resolution_;

  // Cameras for the six cube faces.
  std::array<base::ProjectiveCamera, 6> face_cameras_;

  // Clip-from-world matrices for the six cube faces.
  std::array<ion::math::Matrix4f, 6> clip_from_world_matrices_;

  // Maps of cell indices that are occupied to the entries (position and weight)
  // in that cell. Split into multiple maps so that threads can insert without
  // synchronization.
  std::vector<std::unordered_map<int64, CellEntry>> cell_index_to_entry_;
};

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_BINNING_POINT_CLOUD_BUILDER_H_
