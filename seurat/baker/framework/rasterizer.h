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

#ifndef VR_SEURAT_BAKER_FRAMEWORK_RASTERIZER_H_
#define VR_SEURAT_BAKER_FRAMEWORK_RASTERIZER_H_

#include <memory>
#include <vector>

#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/baker/framework/frame.h"
#include "seurat/baker/framework/ray_bundle.h"
#include "seurat/baker/framework/ray_classifier.h"
#include "seurat/baker/framework/sample_accumulator.h"
#include "seurat/base/status.h"

namespace seurat {
namespace baker {

// Loops over all Views, accumulating their samples onto a set of Frames.
class FrameRasterizer {
 public:
  FrameRasterizer(int thread_count,
                  std::unique_ptr<ingest::ViewGroupLoader> view_loader,
                  std::unique_ptr<RayClassifier> ray_classifier)
      : thread_count_(thread_count),
        view_loader_(std::move(view_loader)),
        ray_classifier_(std::move(ray_classifier)) {}

  // Performs a single pass over all views, accumulating samples.
  base::Status Run(absl::Span<const Frame> frames,
                   absl::Span<const std::shared_ptr<SampleAccumulator>>
                       frame_accumulators) const;

 private:
  // The maximum number of threads to use.
  const int thread_count_;

  // Loads all views of the scene.
  const std::unique_ptr<ingest::ViewGroupLoader> view_loader_;

  // Determines when to use a ray as a solid or freespace sample of a particular
  // Frame.
  const std::unique_ptr<RayClassifier> ray_classifier_;
};

}  // namespace baker
}  // namespace seurat

#endif  // VR_SEURAT_BAKER_FRAMEWORK_RASTERIZER_H_
