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

#include "seurat/artifact/sort_atlas_tiles_transform.h"

#include <algorithm>
#include <vector>

#include "absl/types/span.h"
#include "seurat/geometry/quad_mesh.h"

namespace seurat {
namespace artifact {

using ion::math::Vector2i;

using geometry::IndexedQuad;
using geometry::QuadMesh;
using image::Image4f;

namespace {

void InvertPermutation(absl::Span<const int> permutation,
                       absl::Span<int> inverse_permutation) {
  // Reads in-order read from permutation; scatters to inverse.
  int destination_index = 0;
  for (auto source_index : permutation) {
    inverse_permutation[source_index] = destination_index;
    ++destination_index;
  }
}

// Sorts sizes lexicographically.
bool TileSizeLessThan(const Vector2i& lhs, const Vector2i& rhs) {
  return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]);
}

}  // namespace

base::Status SortAtlasTilesTransform::Process(Artifact* artifact) const {
  CHECK(artifact->quad_mesh);
  geometry::QuadMesh quad_mesh = *artifact->quad_mesh;

  std::vector<int> sort_order(quad_mesh.textures.size());
  std::iota(sort_order.begin(), sort_order.end(), 0);
  std::sort(sort_order.begin(), sort_order.end(),
            [&quad_mesh](const int lhs, const int rhs) {
              return TileSizeLessThan(quad_mesh.textures[lhs].GetSize(),
                                      quad_mesh.textures[rhs].GetSize());
            });

  std::vector<IndexedQuad> sorted_quads(quad_mesh.quads);
  std::vector<Image4f> sorted_textures(quad_mesh.textures.size());

  std::vector<int> inverse_sort_order(sort_order.size());
  InvertPermutation(sort_order, absl::MakeSpan(inverse_sort_order));

  int quad_index = 0;
  for (auto& quad : sorted_quads) {
    // Reorders texture tile images by the calculated sort_order.
    sorted_textures[quad_index] =
        std::move(quad_mesh.textures[sort_order[quad_index]]);
    // Reorders the mapping into the texture set by the inverse sort order.
    // N.B. tex_index(quad_index) === tex_index(p(p^-1(quad_index))).
    quad.texture_index = inverse_sort_order[quad_index];
    ++quad_index;
  }

  artifact->quad_mesh = std::make_shared<QuadMesh>(std::move(sorted_quads),
                                                   std::move(sorted_textures));

  return base::OkStatus();
}

}  // namespace artifact
}  // namespace seurat
