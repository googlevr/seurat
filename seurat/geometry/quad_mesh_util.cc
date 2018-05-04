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

#include "seurat/geometry/quad_mesh_util.h"

#include <array>
#include <vector>

#include "seurat/base/array2d_util.h"
#include "seurat/geometry/cube_face.h"
#include "seurat/geometry/mesh_util.h"
#include "seurat/image/atlaser.h"
#include "seurat/image/image_util.h"

namespace seurat {
namespace geometry {

using geometry::IndexedQuad;
using geometry::Mesh;
using image::Atlaser;
using image::Image4f;
using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2f;
using ion::math::Vector2i;

void MeshFromQuadMesh(absl::Span<const IndexedQuad> rasterized_quads,
                      absl::Span<const Image4f> textures,
                      const Atlaser& atlaser, Mesh* mesh, Image4f* atlas) {
  std::vector<Vector2i> quad_texture_sizes;
  for (const auto& texture_tile : textures) {
    quad_texture_sizes.push_back(texture_tile.GetSize());
  }
  Vector2i total_texture_size;
  std::vector<Point2i> layout(quad_texture_sizes.size());
  atlaser.LayoutTiles(quad_texture_sizes, &total_texture_size,
                      absl::MakeSpan(layout));
  Vector2i atlas_size_target = atlaser.GetAtlasSizeTarget();
  if (0 < atlas_size_target[0]) {
    CHECK_LT(0, atlas_size_target[1]);
    CHECK_LE(total_texture_size[0], atlas_size_target[0]);
    CHECK_LE(total_texture_size[1], atlas_size_target[1]);
    total_texture_size = atlas_size_target;
  }
  // Copy all images into the |atlas| based on the layout computed above.
  if (atlas) {
    atlas->Resize(total_texture_size[0], total_texture_size[1]);
    for (int tile_index = 0; tile_index < textures.size(); ++tile_index) {
      const Vector2i offset = layout[tile_index] - Point2i::Zero();
      base::CopyArray(textures[tile_index], atlas, offset);
    }
  }
  for (int i = 0; i < rasterized_quads.size(); ++i) {
    const IndexedQuad& dq = rasterized_quads[i];
    const Vector2i& tsize = textures[dq.texture_index].GetSize();
    const Vector2i offset = layout[dq.texture_index] - Point2i::Zero();

    // Use half-pixel offsets to compensate for OpenGL's texture
    // parameterization.
    //
    // This ensures that the corners of the quad correspond to pixel centers.
    const std::array<float, 4>& w = dq.texcoord_w;
    std::array<Point3f, 4> texcoords = {
        {{0.5f, 0.5f, 1.0f},
         {(tsize[0] - 0.5f), 0.5f, 1.0f},
         {(tsize[0] - 0.5f), (tsize[1] - 0.5f), 1.0f},
         {0.5f, (tsize[1] - 0.5f), 1.0f}}};
    for (int i = 0; i < 4; ++i) {
      texcoords[i] *= w[i];
    }

    for (auto& pt : texcoords) {
      pt[0] = (pt[0] + offset[0] * pt[2]) / total_texture_size[0];
      pt[1] = (pt[1] + offset[1] * pt[2]) / total_texture_size[1];
    }

    if (mesh) {
      geometry::AppendTriangleFan(dq.quad, {texcoords}, mesh);
    }
  }
}

}  // namespace geometry
}  // namespace seurat
