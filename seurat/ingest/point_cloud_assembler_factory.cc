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

#include "seurat/ingest/point_cloud_assembler_factory.h"

#include "seurat/geometry/binning_point_cloud_builder.h"
#include "seurat/ingest/view_group_loader_util.h"

namespace seurat {
namespace ingest {

using geometry::BinningPointCloudBuilder;
using geometry::PointCloudBuilder;

// static
std::unique_ptr<PointCloudAssembler>
PointCloudAssemblerFactory::CreatePointCloudAssembler(
    std::unique_ptr<ViewGroupLoader> view_group_loader) {
  ion::math::Range3f headbox = ingest::ComputeHeadbox(*view_group_loader);
  ion::math::Vector3i resolution;
  float near_clip;
  BinningPointCloudBuilder::ComputeResolutionAndNearClip(
      headbox, parameters_.pixels_per_degree, &resolution, &near_clip);
  std::unique_ptr<PointCloudBuilder> builder(new BinningPointCloudBuilder(
      parameters_.thread_count, resolution, near_clip));

  return std::unique_ptr<PointCloudAssembler>(new PointCloudAssembler(
      parameters_.thread_count, std::move(view_group_loader),
      std::move(builder)));
}

}  // namespace ingest
}  // namespace seurat
