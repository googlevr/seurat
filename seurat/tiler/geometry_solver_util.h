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

#ifndef VR_SEURAT_TILER_GEOMETRY_SOLVER_UTIL_H_
#define VR_SEURAT_TILER_GEOMETRY_SOLVER_UTIL_H_

#include <utility>
#include <vector>

#include "ceres/ceres.h"
#include "absl/types/span.h"
#include "seurat/tiler/geometry_model.h"
#include "seurat/tiler/point_set.h"

namespace seurat {
namespace tiler {

// Returns a weighted average of the points of interest.
ion::math::Point3f ComputeInitialCenterPoint(
    const PointSet& point_set, absl::Span<const int> points_of_interest);

// Evaluates the |problem| once, returning whether evaluation succeeded.
bool InitialEvaluationSucceeds(ceres::Problem* problem);

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_GEOMETRY_SOLVER_UTIL_H_
