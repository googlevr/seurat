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

#include "seurat/image/image_test_utils.h"

#include "seurat/base/array2d_util.h"
#include "seurat/base/color.h"
#include "seurat/image/image.h"

namespace seurat {
namespace image {

using base::Color4f;
using ion::math::Point2i;
using ion::math::Vector2i;

void StubColorProcessor::ProcessColors(
    absl::Span<base::Color4f> color_set) const {}

Image4f MakeCheckerboard(const Vector2i& size, int square_size,
                         const std::vector<Color4f>& colors) {
  Image4f image(size);
  base::SpatialFillArray(&image, [square_size, &colors](const Point2i& pos) {
    return colors[(pos[0] / square_size + pos[1] / square_size) %
                  colors.size()];
  });
  return image;
}

Image4f ConcatenateImages(const std::vector<Image4f>& images,
                          int separator_width) {
  int height = 0;
  int width = 0;
  for (const auto& image : images) {
    width += image.Width() + separator_width;
    height = std::max(image.Height(), height);
  }
  Image4f result(width, height, Color4f::Zero());
  Vector2i dst_offset(0, 0);
  for (const auto& image : images) {
    base::CopyArray(image, &result, dst_offset);
    dst_offset[0] += image.Width() + separator_width;
  }
  return result;
}

}  // namespace image
}  // namespace seurat
