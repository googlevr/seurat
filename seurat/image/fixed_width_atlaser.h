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

#ifndef VR_SEURAT_IMAGE_FIXED_WIDTH_ATLASER_H_
#define VR_SEURAT_IMAGE_FIXED_WIDTH_ATLASER_H_

#include <vector>

#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/image/atlaser.h"

namespace seurat {
namespace image {

// Implementation of the Atlaser interface that constrains the width of the
// layout to be |texture_width|.
class FixedWidthAtlaser : public Atlaser {
 public:
  explicit FixedWidthAtlaser(const ion::math::Vector2i& atlas_size_target)
      : atlas_size_target_(atlas_size_target) {}

  ion::math::Vector2i GetAtlasSizeTarget() const override {
    return atlas_size_target_;
  }

  void LayoutTiles(absl::Span<const ion::math::Vector2i> tile_sizes,
                   ion::math::Vector2i* total_size,
                   absl::Span<ion::math::Point2i> tile_origins) const override;

 private:
  // Requested size of the texture atlas.
  const ion::math::Vector2i atlas_size_target_;
};

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_FIXED_WIDTH_ATLASER_H_
