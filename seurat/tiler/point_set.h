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

#ifndef VR_SEURAT_TILER_POINT_SET_H_
#define VR_SEURAT_TILER_POINT_SET_H_

#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/base/color.h"

namespace seurat {
namespace tiler {

// An unowning structure-of-arrays of a point cloud, as used by the Tiler.
struct PointSet {
  // A value for PointSet.id which must never be used for a valid PointSet
  // instance.
  static constexpr int kInvalidId = -1;

  // A unique-id specifying this point set.
  //
  // This may be useful for caching expensive preprocessing operations when the
  // same PointSet is being used.
  int id;

  // The positions of all points.
  absl::Span<const ion::math::Point3f> positions;

  // The normals of all points, or empty if no normals are available.
  absl::Span<const ion::math::Vector3f> normals;

  // The colors of all points, or empty if no colors are available.
  absl::Span<const base::Color3f> colors;

  // The weights of all points, or empty if no weights are available.
  //
  // Weights determine how to scale the error-metric evaluated for each point.
  // Higher weight places more importance on accurate reconstruction of those
  // points.
  absl::Span<const float> weights;
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_POINT_SET_H_
