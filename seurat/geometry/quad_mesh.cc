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

#include "seurat/geometry/quad_mesh.h"

namespace seurat {
namespace geometry {

using geometry::Quad3f;
using image::Image4f;

namespace {

bool IndexedQuadsAreValid(const std::vector<IndexedQuad>& quads,
                          const std::vector<Image4f>& textures) {
  int quad_num = 0;
  for (auto const& quad : quads) {
    if (quad.texture_index < 0 || quad.texture_index >= textures.size())
      return false;
    ++quad_num;
  }
  return true;
}

}  // namespace

QuadMesh::QuadMesh(std::vector<IndexedQuad> indexed_quads,
                   std::vector<image::Image4f> indexed_textures)
    : quads(std::move(indexed_quads)), textures(std::move(indexed_textures)) {
  CHECK(IndexedQuadsAreValid(quads, textures))
      << "Quad uses non-existent texture.";
}

}  // namespace geometry
}  // namespace seurat
