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

#ifndef VR_SEURAT_COMPONENT_SHAPE_RENDERABLE_H_
#define VR_SEURAT_COMPONENT_SHAPE_RENDERABLE_H_

#include <memory>

#include "ion/gfx/shape.h"
#include "seurat/component/renderable.h"
#include "seurat/component/shape_component.h"

namespace seurat {
namespace component {

// This class implements the Renderable interface for a ShapeComponent. It
// supports showing arbitrary simple Ion shapes for visualization, debugging,
// and utility rendering.
class ShapeRenderable : public Renderable {
 public:
  static std::unique_ptr<ShapeRenderable> Create(
      const ShapeComponent* shape_component);
  ~ShapeRenderable() override;

  // Renderable implementation.
  Properties GetProperties() const override;
  NodeSet GetRenderNodes(const RenderingContext& context) override;
  void Update(const ViewState& view_state) override;

 private:
  explicit ShapeRenderable(const ShapeComponent* shape_component);

  // Create the shader uniforms.
  void CreateUniforms(const ion::gfx::ShaderInputRegistryPtr& shader_registry);

  // The ShapeComponent which holds the shape state.
  const ShapeComponent* const shape_component_;

  // The Ion node to hold the Shape for drawing the ShapeComponent.
  ion::gfx::NodePtr shape_node_;

  // The Shape which draws this Renderable's set of primitives.
  ion::gfx::ShapePtr shape_;

  // The node set for this renderable.
  NodeSet node_set_;

  // Binds to index of uClipFromShapeMatrix in vertex shader.
  size_t u_leftright_clip_from_shape_matrix_index_;

  DISALLOW_COPY_AND_ASSIGN(ShapeRenderable);
};

}  // namespace component
}  // namespace seurat

#endif  // VR_SEURAT_COMPONENT_SHAPE_RENDERABLE_H_
