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

#ifndef VR_SEURAT_BAKER_FRAMEWORK_SAMPLE_ACCUMULATOR_H_
#define VR_SEURAT_BAKER_FRAMEWORK_SAMPLE_ACCUMULATOR_H_

#include <memory>
#include <vector>

#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/baker/framework/ray_bundle.h"

namespace seurat {
namespace baker {

// Accumulates samples for a single quad.
//
// For example, an implementation may splat color data into a texture, or
// another instance may maintain some representation of the silhouette of the
// Frame.
class SampleAccumulator {
 public:
  virtual ~SampleAccumulator() = default;

  // Adds solid samples & freespace rays.
  //
  // Implementations should merge these into their internal representation, i.e.
  // splat color data, filtering as necessary.
  virtual void Add(
      const RayBundle& bundle,
      absl::Span<const RayBundle::RayIntersectionIndex> solid_samples,
      absl::Span<const int> freespace_rays) = 0;
};

}  // namespace baker
}  // namespace seurat

#endif  // VR_SEURAT_BAKER_FRAMEWORK_SAMPLE_ACCUMULATOR_H_
