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

#include "seurat/baker/framework/rasterizer.h"

#include "seurat/base/parallel.h"
#include "seurat/base/progress.h"
#include "seurat/ingest/view_group_loader_util.h"

namespace seurat {
namespace baker {

using base::Camera;

base::Status FrameRasterizer::Run(
    absl::Span<const Frame> frames,
    absl::Span<const std::shared_ptr<SampleAccumulator>> frame_accumulators)
    const {
  CHECK_EQ(frames.size(), frame_accumulators.size());
  ray_classifier_->Init(frames);

  // Loop over all views from all view-groups.
  base::ScopedProgressRange progress("Generating textures",
                                     view_loader_->GetNumViewGroups());
  return ingest::ForEachViewGroupPrefetching(
      *view_loader_, [&](std::vector<std::shared_ptr<Camera>> cameras,
                         std::vector<image::Ldi4f> ldis) {
        ViewGroupRayBundle bundle(std::move(cameras), std::move(ldis));
        std::vector<RayClassifier::ClassifiedRays> classified_rays =
            ray_classifier_->ClassifyRays(bundle);
        base::ParallelFor(thread_count_, frames.size(), [&](int frame_index) {
          frame_accumulators[frame_index]->Add(
              bundle, classified_rays[frame_index].solid_samples,
              classified_rays[frame_index].freespace_rays);
        });
        progress.IncrementRange(1);
        return base::OkStatus();
      });
}

}  // namespace baker
}  // namespace seurat
