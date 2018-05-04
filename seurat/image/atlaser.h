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

#ifndef VR_SEURAT_IMAGE_ATLASER_H_
#define VR_SEURAT_IMAGE_ATLASER_H_

#include "ion/math/vector.h"
#include "absl/types/span.h"

namespace seurat {
namespace image {

// Interface for organizing multiple smaller textures into a large texture
// atlas.
class Atlaser {
 public:
  virtual ~Atlaser() = default;

  // Returns the target on the atlas size, or (0,0) in case no such target has
  // been set.
  virtual ion::math::Vector2i GetAtlasSizeTarget() const = 0;

  // Lays out a set of tiles, defined by their |tile_sizes|. Returns the
  // |tile_origins| in the texture atlas in a vector. The order of the result
  // vector is the same as the order of the input vector. The final size of the
  // atlas is returned in |total_size|.
  virtual void LayoutTiles(
      absl::Span<const ion::math::Vector2i> tile_sizes,
      ion::math::Vector2i *total_size,
      absl::Span<ion::math::Point2i> tile_origins) const = 0;
};

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_ATLASER_H_
