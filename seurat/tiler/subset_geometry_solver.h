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

#ifndef VR_SEURAT_TILER_SUBSET_GEOMETRY_SOLVER_H_
#define VR_SEURAT_TILER_SUBSET_GEOMETRY_SOLVER_H_

#include "ion/math/vector.h"
#include "seurat/tiler/geometry_model.h"
#include "seurat/tiler/geometry_solver.h"
#include "seurat/tiler/point_set.h"

namespace seurat {
namespace tiler {

// Accelerates any existing GeometrySolver implementation by using only a subset
// of all points for each call to FitModel().
class SubsetGeometrySolver : public GeometrySolver {
 public:
  SubsetGeometrySolver(int max_point_count,
                       std::shared_ptr<GeometrySolver> delegate)
      : max_point_count_(max_point_count), delegate_(std::move(delegate)) {}
  ~SubsetGeometrySolver() override = default;

  // GeometrySolver implementation.
  void Init(const PointSet& point_set) override { delegate_->Init(point_set); }
  void InitializeModel(int point_index, GeometryModel* model) const override {
    delegate_->InitializeModel(point_index, model);
  }
  bool FitModel(absl::Span<const int> point_indices,
                GeometryModel* model) const override;
  float ComputeError(int point_index,
                     const GeometryModel& model) const override {
    return delegate_->ComputeError(point_index, model);
  }

 private:
  const int max_point_count_;

  // The GeometrySolver implementation to accelerate.
  const std::shared_ptr<GeometrySolver> delegate_;
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SUBSET_GEOMETRY_SOLVER_H_
