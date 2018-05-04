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

#ifndef VR_SEURAT_MESH_MESH_RENDERABLE_H_
#define VR_SEURAT_MESH_MESH_RENDERABLE_H_

#include "ion/gfx/shaderinputregistry.h"
#include "seurat/component/renderable.h"
#include "seurat/mesh/mesh_component.h"

namespace seurat {
namespace mesh {

class MeshRenderable : public component::Renderable {
 public:
  // Constructs a mesh renderable.
  explicit MeshRenderable(const MeshComponent* mesh_component);

  Properties GetProperties() const override;
  NodeSet GetRenderNodes(const RenderingContext& context) override;
  void Update(const ViewState& view_state) override;

 private:
  // Pointer to the mesh component.
  const MeshComponent* const mesh_component_;

  // The Ion node for the mesh.
  ion::gfx::NodePtr node_;

  // Uniform index for the array of matrices that transform from component-local
  // space to the viewer's left/right eye clip space.
  size_t u_clip_from_mesh_matrix_index_;
};

}  // namespace mesh
}  // namespace seurat

#endif  // VR_SEURAT_MESH_MESH_RENDERABLE_H_
