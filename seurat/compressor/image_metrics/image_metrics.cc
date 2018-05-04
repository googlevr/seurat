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

#include "seurat/compressor/image_metrics/image_metrics.h"

#include "ion/gfx/image.h"
#include "ion/math/vector.h"

namespace seurat {
namespace compressor {

using ion::math::Vector2i;

float EvalBitrate(const Vector2i& image_size, ion::gfx::Image::Format format) {
  constexpr float kBitsPerByte = 8.0f;
  return kBitsPerByte *
         ion::gfx::Image::ComputeDataSize(format, image_size[0], image_size[1]);
}

}  // namespace compressor
}  // namespace seurat
