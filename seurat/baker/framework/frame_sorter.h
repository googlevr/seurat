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

#ifndef VR_SEURAT_BAKER_FRAMEWORK_FRAME_SORTER_H_
#define VR_SEURAT_BAKER_FRAMEWORK_FRAME_SORTER_H_

#include <functional>
#include <memory>
#include <vector>

#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/baker/framework/frame.h"

namespace seurat {
namespace baker {

// Given a set of Frames, determines their draw order.
//
// Approximate draw order is computed based on the median point of all points
// which would be assigned to a Frame.
class FrameSorter {
 public:
  explicit FrameSorter(int thread_count) : thread_count_(thread_count) {}

  // Overwrites the draw_order of the given |frames| to render back-to-front.
  void ComputeDrawOrder(absl::Span<const ion::math::Point3f> points,
                        absl::Span<Frame> frames);

 private:
  // The maximum number of threads to use.
  const int thread_count_;
};

}  // namespace baker
}  // namespace seurat

#endif  // VR_SEURAT_BAKER_FRAMEWORK_FRAME_SORTER_H_
