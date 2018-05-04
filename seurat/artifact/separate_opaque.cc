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

#include "seurat/artifact/separate_opaque.h"
#include "seurat/image/image_util.h"

namespace seurat {
namespace artifact {

using geometry::IndexedQuad;
using geometry::QuadMesh;
using image::Image4f;

base::Status SeparateOpaque::Process(Artifact* artifact) const {
  CHECK(artifact->quad_mesh)
      << "The input artifact does not support the required representation";
  const QuadMesh& input_quad_mesh = *artifact->quad_mesh;

  std::vector<IndexedQuad> separated_quads;
  std::vector<Image4f> separated_textures;

  // Maps an input texture index to the corresponding output texture index.
  // Needed because having |texture_index| in geometry::IndexedQuad allows for
  // quads and textures appearing in different orders in their respective arrays
  // and also for the same texture being shared by multiple quads.
  std::map<int, int> index_map;

  // Represents the set of textures that have been tested and rejected. When
  // some texture in the |rejected| set is referenced by another quad one can
  // avoid testing it again.
  std::set<int> rejected;

  int output_texture_index = 0;
  for (const auto& quad : input_quad_mesh.quads) {
    int texture_index = quad.texture_index;
    if (0 != rejected.count(texture_index)) continue;
    if (0 == index_map.count(texture_index)) {
      const Image4f& texture = input_quad_mesh.textures[texture_index];
      const bool texture_is_opaque = image::IsOpaque(texture, alpha_threshold_);
      if ((retain_ == Retain::kRetainOpaque) == texture_is_opaque) {
        separated_quads.push_back(
            {quad.quad, quad.texcoord_w, output_texture_index});
        separated_textures.push_back(texture);
        index_map[texture_index] = output_texture_index;
        ++output_texture_index;
      } else {
        rejected.insert(texture_index);
      }
    } else {
      separated_quads.push_back(
          {quad.quad, quad.texcoord_w, index_map[texture_index]});
    }
  }

  artifact->quad_mesh = std::make_shared<QuadMesh>(
      std::move(separated_quads), std::move(separated_textures));

  return base::OkStatus();
}

}  // namespace artifact
}  // namespace seurat
