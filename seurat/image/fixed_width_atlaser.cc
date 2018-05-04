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

#include "seurat/image/fixed_width_atlaser.h"

#include <algorithm>
#include <numeric>

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/base/util.h"
#include "seurat/image/atlaser_util.h"

namespace seurat {
namespace image {

using ion::math::Point2i;
using ion::math::Range2i;
using ion::math::Vector2i;

void FixedWidthAtlaser::LayoutTiles(absl::Span<const Vector2i> tile_sizes,
                                    Vector2i* total_size,
                                    absl::Span<Point2i> tile_origins) const {
  DCHECK_EQ(tile_sizes.size(), tile_origins.size());
  int height;
  ComputeAtlasLayout(tile_sizes, atlas_size_target_[0], &height, tile_origins);
  *total_size = ComputeAtlasSize(tile_sizes, tile_origins);
}

}  // namespace image
}  // namespace seurat
