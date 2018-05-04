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

#ifndef VR_SEURAT_TILER_GEOMETRY_SOLVER_H_
#define VR_SEURAT_TILER_GEOMETRY_SOLVER_H_

#include <vector>

#include "absl/types/span.h"
#include "seurat/tiler/geometry_model.h"
#include "seurat/tiler/point_set.h"

namespace seurat {
namespace tiler {

// An interface for specifying methods for fitting GeometryModels used for
// partitioning.
//
// Both ComputeError() and FitModel() should use the same underlying
// error-metric.  Both convex and non-convex error metrics may be used.
class GeometrySolver {
 public:
  // Performs any initialization required for fitting GeometryModels against the
  // specified PointSet.
  virtual void Init(const PointSet& point_set) = 0;

  // Initializes a GeometryModel which fits the point with the specified index
  // into the PointSet.
  virtual void InitializeModel(int point_index, GeometryModel* model) const = 0;

  // Fits a GeometryModel to the specified points.  The previously-estimated
  // |model| geometry is provided as a starting point, and will be modified by
  // this method.
  //
  // Returns true upon success, false if no model could be fit (e.g. because
  // there are too few points to fit a plane).
  virtual bool FitModel(absl::Span<const int> point_indices,
                        GeometryModel* model) const = 0;

  // Computes an error-metric measuring the deviation of the specified point
  // from the |model|.
  virtual float ComputeError(int point_index,
                             const GeometryModel& model) const = 0;

  virtual ~GeometrySolver() = default;
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_GEOMETRY_SOLVER_H_
