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

#ifndef VR_SEURAT_COMPRESSOR_IMAGE_METRICS_IMAGE_METRICS_H_
#define VR_SEURAT_COMPRESSOR_IMAGE_METRICS_IMAGE_METRICS_H_

#include <cmath>

#include "ion/gfx/image.h"
#include "ion/math/range.h"
#include "seurat/image/image.h"

namespace seurat {
namespace compressor {

// Returns the sum of squared differences between corresponding regions in
// images |a| and |b|. The regions are defined by |range|. The sum of squared
// differences is used as an additive metric for image distortion.
template <typename T>
float EvalRangeSSD(const image::Image<T>& a, const image::Image<T>& b,
                   const ion::math::Range2i& range) {
  CHECK_EQ(a.GetSize(), b.GetSize());
  CHECK_LE(0, range.GetMinPoint()[0]);
  CHECK_LE(0, range.GetMinPoint()[1]);
  CHECK_LT(range.GetMaxPoint()[0], a.GetSize()[0]);
  CHECK_LT(range.GetMaxPoint()[1], a.GetSize()[1]);
  float ssd = 0.0f;
  for (int y = range.GetMinPoint()[1]; y <= range.GetMaxPoint()[1]; ++y) {
    for (int x = range.GetMinPoint()[0]; x <= range.GetMaxPoint()[0]; ++x) {
      const base::Color<T::kDimension, float> color_difference =
          a.At(x, y) - b.At(x, y);
      for (int d = 0; d < T::kDimension; ++d) {
        ssd += color_difference[d] * color_difference[d];
      }
    }
  }
  return ssd;
}

// Returns the sum of squared differences between images |a| and |b|.
template <typename T>
float EvalSSD(const image::Image<T>& a, const image::Image<T>& b) {
  const ion::math::Range2i range = ion::math::Range2i::BuildWithSize(
      ion::math::Point2i::Zero(), a.GetSize() - ion::math::Vector2i(1, 1));
  return EvalRangeSSD(a, b, range);
}

// Returns the peak signal-to-noise ratio corresponding to the difference
// between images |a| and |b|. Both |a| and |b| are assumed to be floating point
// images with channels normalized to the [0.0f, 1.0f] interval.
template <typename T>
float EvalPSNR(const image::Image<T>& a, const image::Image<T>& b) {
  const float ssd = EvalSSD(a, b);
  return 10.0f * std::log10(static_cast<float>(T::kDimension * a.size()) / ssd);
}

// Returns the number of bits used by a four channel RGBA |image| when using
// the Ion image format |format|. |format| may be either kRgba8888 or one of the
// ASTC RGBA formats.
float EvalBitrate(const ion::math::Vector2i& image_size,
                  ion::gfx::Image::Format format);

}  // namespace compressor
}  // namespace seurat

#endif  // VR_SEURAT_COMPRESSOR_IMAGE_METRICS_IMAGE_METRICS_H_
