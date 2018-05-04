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

#include "seurat/ingest/point_cloud_assembler.h"

#include "ion/math/matrix.h"
#include "ion/math/matrixutils.h"
#include "ion/math/transformutils.h"
#include "seurat/base/camera.h"
#include "seurat/base/parallel.h"
#include "seurat/base/progress.h"
#include "seurat/base/util.h"
#include "seurat/image/image.h"
#include "seurat/image/ldi.h"
#include "seurat/ingest/view_group_loader_util.h"

using ion::math::Matrix4f;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2i;
using seurat::base::Camera;
using seurat::image::Image1f;
using seurat::image::Image3f;

namespace seurat {
namespace ingest {

base::Status PointCloudAssembler::Build(std::vector<Point3f>* positions,
                                        std::vector<float>* weights) const {
  ion::math::Range3f view_region;

  const int num_view_groups = view_group_loader_->GetNumViewGroups();
  base::ScopedProgressRange progress("Reading input", num_view_groups);
  SEURAT_RETURN_IF_ERROR(ingest::ForEachViewGroupPrefetching(
      *view_group_loader_, [&](std::vector<std::shared_ptr<Camera>> cameras,
                               std::vector<image::Ldi4f> ldis) {
        for (const auto& camera : cameras) {
          view_region.ExtendByPoint(camera->GetWorldFromEye() *
                                    Point3f::Zero());
        }
        CHECK_EQ(cameras.size(), ldis.size());
        for (int view_index = 0; view_index < cameras.size(); ++view_index) {
          SEURAT_RETURN_IF_ERROR(builder_->AddPoints(PositionsFromView(
              thread_count_, *cameras[view_index], ldis[view_index])));
        }
        progress.IncrementRange(1);
        return base::OkStatus();
      }));

  builder_->GetPositionsAndWeights(positions, weights);
  return base::OkStatus();
}

// static
std::vector<Point3f> PointCloudAssembler::PositionsFromView(
    int thread_count, const base::Camera& camera, const image::Ldi4f& ldi) {
  const Vector2i image_size = ldi.GetSize();
  const int total_sample_count = ldi.GetSampleCount();
  const int pixel_count = image_size[0] * image_size[1];
  std::vector<Point3f> positions(total_sample_count);
  base::ParallelFor(thread_count, pixel_count, [&](int pixel_index) {
    const Point2i pixel(pixel_index % image_size[0],
                        pixel_index / image_size[0]);
    int sample_count = 0;
    int offset = 0;
    ldi.GetSampleCountAndOffset(pixel, &offset, &sample_count);
    for (int sample_index = 0; sample_index < sample_count; ++sample_index) {
      positions[offset + sample_index] =
          camera.RayEnd(pixel, ldi.GetDepths(pixel)[sample_index]);
    }
  });
  return positions;
}

}  // namespace ingest
}  // namespace seurat
