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

#include "seurat/compressor/resampler/box_downsampler.h"

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "seurat/base/array2d_util.h"
#include "seurat/image/image.h"

namespace seurat {
namespace compressor {

using image::Image4f;
using ion::math::Range1f;
using ion::math::Vector2f;
using ion::math::Vector2i;

Image4f BoxDownsampler::Resample(const Image4f& image,
                                 const Vector2i& target_size) const {
  CHECK_NE(target_size, image.GetSize());
  CHECK_LE(2, target_size[0]);
  CHECK_LE(2, target_size[1]);
  CHECK_LE(target_size[0], image.Width());
  CHECK_LE(target_size[1], image.Height());
  Image4f resized_image(target_size);
  // Account for the half pixel border around the texture. The goal is to match
  // as close as possible the regions of the source image and of the downsampled
  // image between texture coordinates [0.5, texture_size-0.5]. This amounts to
  // subtracting 1 from both texture extents when computing the x and y scaling
  // factors.
  Vector2f scaling(
      static_cast<float>(image.Width() - 1) / (target_size[0] - 1),
      static_cast<float>(image.Height() - 1) / (target_size[1] - 1));
  for (int y = 0; y < target_size[1]; ++y) {
    for (int x = 0; x < target_size[0]; ++x) {
      // Range of the box filter within the initial image.
      Range1f x_range((x - 0.5f) * scaling[0], (x + 0.5f) * scaling[0]);
      Range1f y_range((y - 0.5f) * scaling[1], (y + 0.5f) * scaling[1]);

      // Set downsampled image pixel to the average value. The box filter
      // average is clamped to [0,1] and quantized.
      resized_image.At(x, y) = ComputeBoxAverage(image, x_range, y_range)
                                   .ClampToUnit()
                                   .AsColorUI8()
                                   .AsColorF();
    }
  }
  return resized_image;
}

}  // namespace compressor
}  // namespace seurat
