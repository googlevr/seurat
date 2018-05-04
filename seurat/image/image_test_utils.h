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

#ifndef VR_SEURAT_IMAGE_IMAGE_TEST_UTILS_H_
#define VR_SEURAT_IMAGE_IMAGE_TEST_UTILS_H_

#include <string>
#include <vector>

#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/image/atlaser.h"
#include "seurat/image/color_processor.h"
#include "seurat/image/image.h"

namespace seurat {
namespace image {

// This is a ColorProcessor implementation that does nothing.
class StubColorProcessor : public ColorProcessor {
 public:
  StubColorProcessor() = default;
  ~StubColorProcessor() override = default;

  // ColorProcessor implementation.
  void ProcessColors(absl::Span<base::Color4f> color_set) const override;
};

// Do-nothing implementation of the Atlaser interface.
class FakeAtlaser : public image::Atlaser {
 public:
  ion::math::Vector2i GetAtlasSizeTarget() const override {
    return ion::math::Vector2i::Zero();
  }

  void LayoutTiles(absl::Span<const ion::math::Vector2i> tile_sizes,
                   ion::math::Vector2i* total_size,
                   absl::Span<ion::math::Point2i> tile_origins) const override {
    *total_size = ion::math::Vector2i::Zero();
  }
};

// Returns an image of size |size| that has the pattern of a checkerboard with
// squares of size |square_size|.
Image4f MakeCheckerboard(const ion::math::Vector2i& size, int square_size,
                         const std::vector<base::Color4f>& colors);

// Returns an image built by concatenating an array of |images|.
Image4f ConcatenateImages(const std::vector<Image4f>& images,
                          int separator_width = 0);

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_IMAGE_TEST_UTILS_H_
