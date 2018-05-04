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

#ifndef VR_SEURAT_GEOMETRY_POINT_CLOUD_BUILDER_H_
#define VR_SEURAT_GEOMETRY_POINT_CLOUD_BUILDER_H_

#include <vector>

#include "ion/math/vector.h"
#include "seurat/base/status.h"

namespace seurat {
namespace geometry {

class PointCloudBuilder {
 public:
  virtual ~PointCloudBuilder() = default;

  // Merges the |positions| into the point cloud.
  virtual base::Status AddPoints(
      const std::vector<ion::math::Point3f>& positions) = 0;

  // Merges the |positions| together with the corresponding point |weights| into
  // the point cloud.
  virtual base::Status AddPointsWithWeights(
      const std::vector<ion::math::Point3f>& positions,
      const std::vector<float>& weights) = 0;

  // Creates a vector with the world-space positions of all points accumulated
  // and a vector with the corresponding weights. Call this method once and only
  // once when all views were added.
  virtual void GetPositionsAndWeights(
      std::vector<ion::math::Point3f>* positions,
      std::vector<float>* weights) = 0;
};

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_POINT_CLOUD_BUILDER_H_
