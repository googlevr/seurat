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

#ifndef VR_SEURAT_IMAGE_NEARLY_SQUARE_ATLASER_H_
#define VR_SEURAT_IMAGE_NEARLY_SQUARE_ATLASER_H_

#include <vector>

#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/image/atlaser.h"

namespace seurat {
namespace image {

// Implementation of the Atlaser interface that produces a texture atlas that is
// approximately a square.
class NearlySquareAtlaser : public Atlaser {
 public:
  ion::math::Vector2i GetAtlasSizeTarget() const override {
    return ion::math::Vector2i::Zero();
  }

  void LayoutTiles(absl::Span<const ion::math::Vector2i> tile_sizes,
                   ion::math::Vector2i* total_size,
                   absl::Span<ion::math::Point2i> tile_origins) const override;
};

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_NEARLY_SQUARE_ATLASER_H_
