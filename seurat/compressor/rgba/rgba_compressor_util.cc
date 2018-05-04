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

#include "seurat/compressor/rgba/rgba_compressor_util.h"

#include <memory>

#include "ion/math/vector.h"
#include "seurat/compressor/resampler/box_downsampler.h"
#include "seurat/compressor/resampler/gl_linear_upsampler.h"
#include "seurat/compressor/rgba/rgba_compressor.h"
#include "seurat/compressor/rgba/rgba_encoding_selector.h"
#include "seurat/compressor/rgba/rgba_rate_resizer.h"
#include "seurat/compressor/rgba/rgba_selecting_compressor.h"
#include "seurat/image/atlaser.h"

namespace seurat {
namespace compressor {

using image::Atlaser;
using ion::math::Vector2i;

std::unique_ptr<RgbaCompressor> BuildAtlasSizeTargetRgbaCompressor(
    int thread_count, const Vector2i& block_size,
    std::shared_ptr<Atlaser> atlaser) {
  Vector2i atlas_size_target = atlaser->GetAtlasSizeTarget();
  CHECK_LT(0, atlas_size_target[0]);
  CHECK_LT(0, atlas_size_target[1]);
  auto build_upsampler = [&]() {
    return std::unique_ptr<Resampler>(new GlLinearUpsampler());
  };
  auto build_downsampler = [&]() {
    return std::unique_ptr<Resampler>(new BoxDownsampler());
  };
  auto build_resizer = [&]() {
    return std::unique_ptr<TextureEncoder>(new RgbaRateResizer(
        thread_count, build_upsampler(), build_downsampler(), block_size));
  };
  auto build_encoding_selector = [&]() {
    return std::unique_ptr<RgbaEncodingSelector>(
        new AtlasSizeTargetSelector(atlaser));
  };
  return std::unique_ptr<RgbaCompressor>(
      new RgbaSelectingCompressor(build_resizer(), build_encoding_selector()));
}

}  // namespace compressor
}  // namespace seurat
