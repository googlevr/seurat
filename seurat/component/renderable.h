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

#ifndef VR_SEURAT_COMPONENT_RENDERABLE_H_
#define VR_SEURAT_COMPONENT_RENDERABLE_H_

#include <vector>

#include "ion/gfx/node.h"
#include "ion/gfx/renderer.h"
#include "ion/gfx/shaderinputregistry.h"
#include "ion/gfxutils/shadermanager.h"
#include "ion/math/matrix.h"
#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "ion/remote/remoteserver.h"

namespace seurat {
namespace component {

// An interface for rendering a Seurat Component.
//
// Renderables are retrieved from Component instances.  Rendering may be a
// multi-pass operation, thus the Renderable exports an ion::gfx::Node instance
// for each pass required, which the Ion renderer can render in order.
class Renderable {
 public:
  // Static properties of the Renderable.
  struct Properties {
    // The axis-aligned bounding box for geometry rendered by this Renderable,
    // in the Renderable's local coordinate space.
    ion::math::Range3f bounding_box;
  };

  struct RenderingContext {
    ion::gfx::RendererPtr renderer;
    ion::gfxutils::ShaderManagerPtr shader_manager;
    ion::math::Vector2i render_target_size;
  };

  // The current state of the VR viewer.
  struct ViewState {
    // Transforms from the component-space to the clip-space of the left eye.
    // The root component is in the start-space.
    ion::math::Matrix4f left_clip_from_component_matrix;
    // Transforms from the component-space to the clip-space of the right eye.
    ion::math::Matrix4f right_clip_from_component_matrix;
    // Transforms from the component-space to the eye-space of the left eye.
    ion::math::Matrix4f left_eye_from_component_matrix;
    // Transforms from the component-space to the eye-space of the right eye.
    ion::math::Matrix4f right_eye_from_component_matrix;
  };

  // Contains (optional) ion nodes to plug into each render pass.
  struct NodeSet {
    // Nodes used by the SingleResolveRenderer
    ion::gfx::NodePtr opaque;
    ion::gfx::NodePtr transparent;

    // Returns pointers to all nodes, in order.
    std::vector<ion::gfx::NodePtr*> AllNodes() {
      return std::vector<ion::gfx::NodePtr*>{&opaque, &transparent};
    }
  };

 public:
  // Returns static properties.
  virtual Properties GetProperties() const = 0;

  // Returns a set of ion nodes to render this component.
  virtual NodeSet GetRenderNodes(const RenderingContext& context) = 0;

  // Updates the nodes returned by |GetRenderNodes| with the given |view_state|.
  virtual void Update(const ViewState& view_state) = 0;

  virtual ~Renderable() = default;
};

}  // namespace component
}  // namespace seurat

#endif  // VR_SEURAT_COMPONENT_RENDERABLE_H_
