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

#ifndef VR_SEURAT_BAKER_FRAMEWORK_SILHOUETTE_ACCUMULATOR_H_
#define VR_SEURAT_BAKER_FRAMEWORK_SILHOUETTE_ACCUMULATOR_H_

#include <memory>

#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/baker/framework/frame.h"
#include "seurat/baker/framework/ray_bundle.h"
#include "seurat/baker/framework/sample_accumulator.h"
#include "seurat/baker/framework/silhouette.h"
#include "seurat/image/image.h"

namespace seurat {
namespace baker {

// Adaptor class to wrap a SilhouetteBuffer as a SampleAccumulator for use with
// the FrameRasterizer.
class SilhouetteAccumulator : public SampleAccumulator {
 public:
  SilhouetteAccumulator(const Frame* frame, int supersample_factor,
                        std::unique_ptr<SilhouetteBuffer> silhouette_buffer)
      : frame_(frame),
        supersample_factor_(supersample_factor),
        silhouette_buffer_(std::move(silhouette_buffer)) {}
  ~SilhouetteAccumulator() override = default;

  // Resolves the silhouette to the |image| at its given resolution.
  void Resolve(image::Image1f* alpha) const;

  // SampleAccumulator implementation.
  void Add(const RayBundle& bundle,
           absl::Span<const std::tuple<int, int>> solid_samples,
           absl::Span<const int> freespace_rays) override;

 private:
  const Frame* frame_;

  // The amount to supersample in both directions.
  //
  // In other words, the total number of samples per pixel is
  // supersample_factor_^2.
  const int supersample_factor_;

  // Represents the 2D silhouette.
  const std::unique_ptr<SilhouetteBuffer> silhouette_buffer_;
};

}  // namespace baker
}  // namespace seurat

#endif  // VR_SEURAT_BAKER_FRAMEWORK_SILHOUETTE_ACCUMULATOR_H_
