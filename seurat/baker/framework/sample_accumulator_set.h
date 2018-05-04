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

#ifndef VR_SEURAT_BAKER_FRAMEWORK_SAMPLE_ACCUMULATOR_SET_H_
#define VR_SEURAT_BAKER_FRAMEWORK_SAMPLE_ACCUMULATOR_SET_H_

#include <memory>
#include <vector>

#include "seurat/baker/framework/rasterizer.h"

namespace seurat {
namespace baker {

// Wraps a set of SampleAccumulators as a single accumulator.
class SampleAccumulatorSet : public SampleAccumulator {
 public:
  explicit SampleAccumulatorSet(
      absl::Span<const std::shared_ptr<SampleAccumulator>> accumulators)
      : accumulators_(accumulators.begin(), accumulators.end()) {}
  ~SampleAccumulatorSet() override = default;

  // SampleAccumulator implementation.
  void Add(const RayBundle& bundle,
           absl::Span<const std::tuple<int, int>> solid_samples,
           absl::Span<const int> freespace_rays) override {
    for (const auto& delegate : accumulators_) {
      delegate->Add(bundle, solid_samples, freespace_rays);
    }
  }

 private:
  std::vector<std::shared_ptr<SampleAccumulator>> accumulators_;
};

}  // namespace baker
}  // namespace seurat

#endif  // VR_SEURAT_BAKER_FRAMEWORK_SAMPLE_ACCUMULATOR_SET_H_
