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

#ifndef VR_SEURAT_INGEST_POINT_CLOUD_ASSEMBLER_H_
#define VR_SEURAT_INGEST_POINT_CLOUD_ASSEMBLER_H_

#include <memory>

#include "ion/math/vector.h"
#include "json/json.h"
#include "seurat/geometry/point_cloud_builder.h"
#include "seurat/image/ldi.h"
#include "seurat/ingest/view_group_loader.h"

namespace seurat {
namespace ingest {

// PointCloudAssembler converts structured view data into a point cloud in
// parallel.
class PointCloudAssembler {
 public:
  PointCloudAssembler(int thread_count,
                      std::unique_ptr<ViewGroupLoader> view_group_loader,
                      std::unique_ptr<geometry::PointCloudBuilder> builder)
      : thread_count_(thread_count),
        view_group_loader_(std::move(view_group_loader)),
        builder_(std::move(builder)) {}

  // Creates the vector of world-space positions.
  base::Status Build(std::vector<ion::math::Point3f>* positions,
                     std::vector<float>* weights) const;

  // Computes and returns a vector with the world-space positions of all points
  // in the view group represented by |cameras| and |ldis|.
  static std::vector<ion::math::Point3f> PositionsFromView(
      int thread_count, const base::Camera& camera, const image::Ldi4f& ldi);

 private:
  // Maximum number of threads to use.
  int thread_count_;

  // The source of color, depth, and camera data.
  std::unique_ptr<ViewGroupLoader> view_group_loader_;

  // The point cloud builder used for the conversion.
  std::unique_ptr<geometry::PointCloudBuilder> builder_;
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_POINT_CLOUD_ASSEMBLER_H_
