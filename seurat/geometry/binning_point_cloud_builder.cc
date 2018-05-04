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

#include "seurat/geometry/binning_point_cloud_builder.h"

#include <algorithm>
#include <limits>
#include <random>

#include "ion/math/matrix.h"
#include "ion/math/matrixutils.h"
#include "ion/math/transformutils.h"
#include "seurat/base/color.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/base/parallel.h"
#include "seurat/base/reporting.h"
#include "seurat/base/util.h"
#include "seurat/geometry/cube_face.h"

namespace seurat {
namespace geometry {

using base::Color3f;
using geometry::CubeFace;
using ion::math::Matrix4f;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Point3i;
using ion::math::Point4f;
using ion::math::Vector2i;
using ion::math::Vector3i;

BinningPointCloudBuilder::BinningPointCloudBuilder(int thread_count,
                                                   const Vector3i& resolution,
                                                   float near_clip)
    : thread_count_(thread_count),
      resolution_(resolution),
      cell_index_to_entry_(thread_count) {
  // Set up the cameras and clip-from-world matrices as six cube-face cameras at
  // the origin (reference cameras). The clip_from_eye matrix uses a
  // far-clipping-plane at infinity.
  const Matrix4f clip_from_eye(              //
      1.0f, 0.0f, 0.0f, 0.0f,                //
      0.0f, 1.0f, 0.0f, 0.0f,                //
      0.0f, 0.0f, -1.0f, -2.0f * near_clip,  //
      0.0f, 0.0f, -1.0f, 0.0f);
  for (int cube_face_index = 0; cube_face_index < 6; ++cube_face_index) {
    const CubeFace cube_face = static_cast<CubeFace>(cube_face_index);
    face_cameras_[cube_face_index] =
        base::ProjectiveCamera(Vector2i(resolution_[0], resolution_[1]),
                               clip_from_eye, LookAtMatrixFromFace(cube_face));
    clip_from_world_matrices_[cube_face_index] =
        face_cameras_[cube_face_index].GetClipFromEye() *
        face_cameras_[cube_face_index].GetEyeFromWorld();
  }
}

void BinningPointCloudBuilder::ComputeResolutionAndNearClip(
    const ion::math::Range3f& headbox, float pixels_per_degree,
    ion::math::Vector3i* resolution, float* near_clip) {
  // Set near clipping plane to distance from headbox center to closest headbox
  // face.
  *near_clip = std::numeric_limits<float>::max();
  const ion::math::Vector3f headbox_size = headbox.GetSize();
  for (int d = 0; d < 3; ++d) {
    if (headbox_size[d] * 0.5f < *near_clip) {
      *near_clip = headbox_size[d] * 0.5f;
    }
  }
  DCHECK(*near_clip < std::numeric_limits<float>::max());
  if (*near_clip == 0.0f) {
    base::SeuratError(
        "Headbox has zero extent along at least one dimension. Arbitrarily "
        "setting binning near clipping plane to 0.1 units from headbox "
        "center.");
    *near_clip = 0.1f;
  }

  // The Z resolution is set to half the X and Y resolution. The Z-range of the
  // last bin extends from -N*Rz to infinity, where N is the near-clip value and
  // Rz is the Z-resolution of the binning grid. If we set N=W, where W is the
  // radius of the headbox, then the last bin covers exactly one pixel of a
  // camera at (0, 0, W) with image resolution Rx, if Rz = Rx/2.
  const int kResolutionScale = 2.0f;
  const int pixels_per_90_degree =
      static_cast<int>(90.0f * pixels_per_degree * kResolutionScale);
  *resolution = ion::math::Vector3i(pixels_per_90_degree, pixels_per_90_degree,
                                    pixels_per_90_degree / 2);
}

base::Status BinningPointCloudBuilder::CellIndexFromPoint(
    const Point3f& position_world, int64* cell_index) const {
  // Find the cube face for the point and transform it into clip space.
  const CubeFace cube_face = geometry::CubeFaceFromPosition(position_world);
  // This must be int64, otherwise the cell_index computation below will be
  // carried out incorrectly in 32-bit, because the resolution vector is only
  // 32-bit.
  const int64 cube_face_index = static_cast<int>(cube_face);
  const ion::math::Matrix4f& clip_from_world =
      clip_from_world_matrices_[cube_face_index];
  const Point4f position_clip_hom =
      clip_from_world * base::Point4FromPoint3(position_world);
  // Check if the point is behind the cube-face camera.
  if (position_clip_hom[3] <= 0.0f) {
    return base::InvalidArgumentError(
        "Point is behind the camera used for binning.  Possible causes for "
        "this problem: geometry inside the headbox, incorrect matrices, "
        "incorrect depth values, other errors in the generation process. "
        "(Showing only the first error)");
  }
  // Check if the point is outside the frustum that defines the perspective grid
  // for binning.
  if (position_clip_hom[0] < -position_clip_hom[3] ||
      position_clip_hom[0] > position_clip_hom[3] ||
      position_clip_hom[1] < -position_clip_hom[3] ||
      position_clip_hom[1] > position_clip_hom[3] ||
      position_clip_hom[2] < -position_clip_hom[3] ||
      position_clip_hom[2] > position_clip_hom[3]) {
    return base::InvalidArgumentError(
        "Point is outside of the frustum used for binning.  Possible causes "
        "for this problem: geometry inside the headbox, incorrect matrices, "
        "incorrect depth values, other errors in the generation process. "
        "(Showing only the first error)");
  }
  const Point3f position_clip = base::Point3FromPoint4(position_clip_hom);

  // Map the position to the discrete perspective grid.
  const ion::math::Point<3, int64> position_clip_discrete(
      (position_clip[0] + 1.0f) * 0.5f * resolution_[0],
      (position_clip[1] + 1.0f) * 0.5f * resolution_[1],
      (position_clip[2] + 1.0f) * 0.5f * resolution_[2]);

  // Compute an index that is unique across the six perspective grids. The
  // computation must be carried out in 64-bit.
  const int64 index =
      position_clip_discrete[0] +                                    //
      position_clip_discrete[1] * resolution_[0] +                   //
      position_clip_discrete[2] * resolution_[0] * resolution_[1] +  //
      cube_face_index * resolution_[0] * resolution_[1] * resolution_[2];

  if (index < 0) {
    return base::OutOfRangeError(
        "Negative binning cell index.  This should never happen.");
  }

  *cell_index = index;
  return base::OkStatus();
}

base::Status BinningPointCloudBuilder::AddPoints(
    const std::vector<Point3f>& positions) {
  constexpr float kDefaultWeight = 1.0f;
  const std::vector<float> weights(positions.size(), kDefaultWeight);
  return AddPointsWithWeights(positions, weights);
}

base::Status BinningPointCloudBuilder::AddPointsWithWeights(
    const std::vector<Point3f>& positions, const std::vector<float>& weights) {
  // The algorithm below executes in two stages. The first stage parallelizes
  // over all input positions. It identifies the points that need to be inserted
  // into the point cloud. The cell indices and entries (positions, weights,
  // etc.) are collected in a two dimensional array of deques. Each thread owns
  // a vector of deques, one for each map of cell indices. In the second stage,
  // we parallelize over the maps of cell indices and insert all points and
  // entries that were collected for this map. If multiple points map to the
  // same cell, we randomly pick one of them. Using reservoir sampling, this
  // works across calls to AddPoints.

  std::vector<std::vector<std::deque<std::pair<int64, CellEntry>>>>
      cell_indices_to_insert(thread_count_);
  for (int tid = 0; tid < thread_count_; ++tid) {
    cell_indices_to_insert[tid].resize(thread_count_);
  }

  std::vector<base::Status> per_thread_status(thread_count_);
  base::ParallelFor(thread_count_, thread_count_, [&](int tid) {
    const int begin = tid * (positions.size() / thread_count_);
    const int end = std::min(positions.size(),
                             (tid + 1) * (positions.size() / thread_count_));
    for (int i = begin; i < end; ++i) {
      int64 cell_index;
      per_thread_status[tid] = CellIndexFromPoint(positions[i], &cell_index);
      if (!per_thread_status[tid]) {
        // Early out on this thread, if an error occurred.
        return;
      }
      // The index into cell_index_to_entry_ where this cell_index would be
      // stored.
      const int map_index = cell_index % thread_count_;
      cell_indices_to_insert[tid][map_index].push_back(
          std::make_pair(cell_index, CellEntry(positions[i], weights[i])));
    }
  });
  // If errors occurred, return the first error from the first thread that
  // failed.
  for (int tid = 0; tid < thread_count_; ++tid) {
    SEURAT_RETURN_IF_ERROR(per_thread_status[tid]);
  }

  // Parallelize over the maps of cell indices.
  base::ParallelFor(thread_count_, thread_count_, [&](int map_index) {
    std::mt19937 random;
    random.seed(map_index);

    for (int tid = 0; tid < thread_count_; ++tid) {
      for (const auto& cell_index_and_entry :
           cell_indices_to_insert[tid][map_index]) {
        auto iter =
            cell_index_to_entry_[map_index].find(cell_index_and_entry.first);
        // Use reservoir sampling to randomly select one of all the points
        // that map to this cell. See Chao, M.T. "A general purpose unequal
        // probability sampling plan", Biometrika (1982), vol. 69, no. 3, pp.
        // 653-6.
        if (iter != cell_index_to_entry_[map_index].end()) {
          iter->second.weight += cell_index_and_entry.second.weight;
          std::bernoulli_distribution dist(cell_index_and_entry.second.weight /
                                           iter->second.weight);
          if (dist(random)) {
            iter->second.position = cell_index_and_entry.second.position;
          }
        } else {
          cell_index_to_entry_[map_index].insert(cell_index_and_entry);
        }
      }
    }
  });

  return base::OkStatus();
}

void BinningPointCloudBuilder::GetPositionsAndWeights(
    std::vector<ion::math::Point3f>* positions, std::vector<float>* weights) {
  CHECK_EQ(thread_count_, cell_index_to_entry_.size());
  int64 num_points = 0;
  std::vector<int> offsets(thread_count_);
  for (int tid = 0; tid < thread_count_; ++tid) {
    offsets[tid] = num_points;
    num_points += cell_index_to_entry_[tid].size();
  }

  positions->clear();
  positions->resize(num_points);
  weights->clear();
  weights->resize(num_points);

  base::ParallelFor(thread_count_, thread_count_, [&](int tid) {
    int offset = offsets[tid];
    for (const auto& cell_index_and_entry : cell_index_to_entry_[tid]) {
      const Point3f& point = cell_index_and_entry.second.position;
      positions->at(offset) = point;
      weights->at(offset) = cell_index_and_entry.second.weight;
      ++offset;
    }
  });

  for (int tid = 0; tid < cell_index_to_entry_.size(); ++tid) {
    cell_index_to_entry_[tid].clear();
  }
}

}  // namespace geometry
}  // namespace seurat
