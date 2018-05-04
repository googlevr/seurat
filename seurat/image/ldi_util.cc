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

#include "seurat/image/ldi_util.h"

#include "ion/math/vector.h"
#include "seurat/base/color.h"

namespace seurat {
namespace image {

using base::Color3f;
using base::Color4f;
using ion::math::Point2i;

Image3f FlattenLdi(const Ldi4f& ldi) {
  Image3f flat_image(ldi.GetSize());
  for (int y = 0; y < ldi.GetHeight(); ++y) {
    for (int x = 0; x < ldi.GetWidth(); ++x) {
      Color3f dst = Color3f::Zero();
      const Point2i pixel(x, y);
      for (int s = ldi.GetSampleCount(pixel) - 1; s >= 0; --s) {
        const Color4f& src = ldi.GetColors(pixel)[s];
        const Color3f src_rgb(src[0], src[1], src[2]);
        const float src_alpha = src[3];
        dst = (1.0f - src_alpha) * dst + src_rgb;
      }
      flat_image.At(pixel) = dst;
    }
  }
  return flat_image;
}

Ldi4f FlipLdiVertically(const Ldi4f& src) {
  const int pixel_count = src.GetWidth() * src.GetHeight();
  std::vector<int> dst_sample_counts;
  std::vector<Color4f> dst_colors;
  std::vector<float> dst_depths;
  dst_sample_counts.reserve(pixel_count);
  dst_colors.reserve(src.GetSampleCount());
  dst_depths.reserve(src.GetSampleCount());

  for (int y = src.GetHeight() - 1; y >= 0; --y) {
    for (int x = 0; x < src.GetWidth(); ++x) {
      const ion::math::Point2i src_pixel(x, y);
      dst_sample_counts.push_back(src.GetSampleCount(src_pixel));
      for (int s = 0; s < src.GetSampleCount(src_pixel); ++s) {
        dst_colors.push_back(src.GetColors(src_pixel)[s]);
        dst_depths.push_back(src.GetDepths(src_pixel)[s]);
      }
    }
  }
  return Ldi4f(src.GetSize(), dst_sample_counts, std::move(dst_colors),
               std::move(dst_depths));
}

}  // namespace image
}  // namespace seurat
