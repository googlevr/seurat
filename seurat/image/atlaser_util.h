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

#ifndef VR_SEURAT_IMAGE_ATLASER_UTIL_H_
#define VR_SEURAT_IMAGE_ATLASER_UTIL_H_

#include "ion/math/vector.h"
#include "absl/types/span.h"

namespace seurat {
namespace image {

// Computes the layout of a set of tiles, defined by their |tile_sizes|. Returns
// the origins of the tiles in the texture atlas in a vector. The order of the
// result vector is the same as the order of the input vector. The layout is
// performed with the given |atlas_width|. The required |atlas_height| is
// returned in a pointer.
void ComputeAtlasLayout(absl::Span<const ion::math::Vector2i> tile_sizes,
                        int atlas_width, int* height,
                        absl::Span<ion::math::Point2i> tile_origins);

// Returns the atlas size given |tile_sizes| and |tile_origins|.
ion::math::Vector2i ComputeAtlasSize(
    absl::Span<const ion::math::Vector2i> tile_sizes,
    absl::Span<const ion::math::Point2i> tile_origins);
}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_ATLASER_UTIL_H_
