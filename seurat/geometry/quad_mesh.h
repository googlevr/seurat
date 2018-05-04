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

#ifndef VR_SEURAT_GEOMETRY_QUAD_MESH_H_
#define VR_SEURAT_GEOMETRY_QUAD_MESH_H_

#include <vector>

#include "seurat/geometry/quad.h"
#include "seurat/image/image.h"

namespace seurat {
namespace geometry {

// Defines a textured (via an index into the owning QuadMesh's texture tile set)
// quadrilateral for one face of a QuadMesh. Intended for use only with
// QuadMesh.
struct IndexedQuad {
  // Constructs an index quad with all values default initialized.
  IndexedQuad() = default;

  // Constructs an indexed quad from all members.
  IndexedQuad(const Quad3f& quad, const std::array<float, 4>& texcoord_w,
              int texture_index)
      : quad(quad), texcoord_w(texcoord_w), texture_index(texture_index) {}

  // The quad is textured such that quad-space coordinates (0, 0) correspond
  // to the center of the texture sample with value texture.At(0, 0).
  geometry::Quad3f quad;

  // The w-component of the texture coordinates for each quad vertex.
  std::array<float, 4> texcoord_w;

  // Index of texture tile in the texture vector, i.e. QuadMesh::textures.
  int texture_index;

  bool operator==(const IndexedQuad& rhs) const {
    return quad == rhs.quad && texcoord_w == rhs.texcoord_w &&
           texture_index == rhs.texture_index;
  }
  bool operator!=(const IndexedQuad& rhs) const { return !operator==(rhs); }
};

// Defines a textured mesh via a set of independent, textured quads.
struct QuadMesh {
  // Constructs an empty quad mesh.
  QuadMesh() = default;

  // Constructs a quad mesh from quads textured by indexing a texture set.
  QuadMesh(std::vector<IndexedQuad> indexed_quads,
           std::vector<image::Image4f> indexed_textures);

  std::vector<IndexedQuad> quads;
  std::vector<image::Image4f> textures;
};

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_QUAD_MESH_H_
