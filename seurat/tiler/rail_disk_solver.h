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

#ifndef VR_SEURAT_TILER_RAIL_DISK_SOLVER_H_
#define VR_SEURAT_TILER_RAIL_DISK_SOLVER_H_

#include "ion/math/vector.h"
#include "seurat/tiler/geometry_model.h"
#include "seurat/tiler/geometry_solver.h"
#include "seurat/tiler/point_set.h"
#include "seurat/tiler/subdivision.h"

namespace seurat {
namespace tiler {

// A solver which treats each GeometryModel instance as a disk, with two
// components to the cost function:
//  1. A "tangential term" penalizing the distance from the center of the disk
//  to the projection of points onto the disk's plane.
//  2. A "normal term" penalizing the distance from a point to its projection on
//  the disk.
//
// In both cases, the projection is performed by intersecting the ray from the
// origin to the point with the disk's plane.
//
// In addition, a soft-constraint is added that the plane of the disk intersect
// the rails of the Subdivision cell in which it lives.
class RailDiskSolver : public GeometrySolver {
 public:
  explicit RailDiskSolver(float tangential_factor,
                          std::shared_ptr<Subdivision> subdivision,
                          ion::math::Range1f depth_range)
      : tangential_factor_(tangential_factor),
        subdivision_(std::move(subdivision)),
        depth_range_(depth_range) {}
  ~RailDiskSolver() override = default;

  // GeometrySolver implementation.
  void Init(const PointSet& point_set) override;
  void InitializeModel(int point_index, GeometryModel* model) const override;
  bool FitModel(absl::Span<const int> point_indices,
                GeometryModel* model) const override;
  float ComputeError(int point_index,
                     const GeometryModel& model) const override;

 private:
  // Scales the tangential term of the cost function.
  const float tangential_factor_;

  // Used to enforce valid rail depths.
  const std::shared_ptr<Subdivision> subdivision_;

  // The range of "normal" depths.  Rail depths beyond this are penalized.
  const ion::math::Range1f depth_range_;

  // The current PointSet in use.
  PointSet point_set_;
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_RAIL_DISK_SOLVER_H_
