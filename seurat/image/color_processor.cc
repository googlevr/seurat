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

#include "seurat/image/color_processor.h"

#include <cmath>

namespace seurat {
namespace image {

void GammaToneMapper::ProcessColors(absl::Span<base::Color4f> color_set) const {
  const float reciprocal_gamma = 1.0f / gamma_;
  for (auto& color : color_set) {
    base::Color4f mapped_color = color;
    if (reciprocal_gamma != 1.0f) {
      mapped_color = base::Color4f(std::pow(mapped_color[0], reciprocal_gamma),
                                   std::pow(mapped_color[1], reciprocal_gamma),
                                   std::pow(mapped_color[2], reciprocal_gamma),
                                   mapped_color[3]);
    }
    color = mapped_color;
  }
}

void PremultipliedAlphaConverter::ProcessColors(
    absl::Span<base::Color4f> color_set) const {
  for (auto& color : color_set) {
    const float alpha = color[3];
    color[0] *= alpha;
    color[1] *= alpha;
    color[2] *= alpha;
  }
}

void ColorProcessorPipeline::ProcessColors(
    absl::Span<base::Color4f> color_set) const {
  for (const auto& processor : sequence_) {
    processor->ProcessColors(color_set);
  }
}

}  // namespace image
}  // namespace seurat
