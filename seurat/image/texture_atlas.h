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

#ifndef VR_SEURAT_IMAGE_TEXTURE_ATLAS_H_
#define VR_SEURAT_IMAGE_TEXTURE_ATLAS_H_

#include <vector>

#include "ion/math/vector.h"
#include "absl/types/span.h"

namespace seurat {
namespace image {

class TextureAtlas {
 public:
  // Calculates and returns a partition of the tiles in |tile_sizes| into
  // contiguous subsets which fit into layouts no larger than |max_size|. The
  // returned partition indicates indices to split |tile_sizes|. The partitions
  // exclude the tile indicated by each split index, so if the resulting vector
  // is {10, 35, 72}, then the first partition is [0, 10), the second [10, 35),
  // and so on.
  static std::vector<int> FindAtlasSplits(
      absl::Span<const ion::math::Vector2i> tile_sizes,
      const ion::math::Vector2i& max_size);
};

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_TEXTURE_ATLAS_H_
