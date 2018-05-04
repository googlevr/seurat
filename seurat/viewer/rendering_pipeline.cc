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

#include "seurat/viewer/rendering_pipeline.h"

#include "ion/base/logging.h"
#include "ion/gfx/graphicsmanager.h"
#include "ion/gfx/renderer.h"

namespace seurat {
namespace viewer {

using component::Component;
using component::Renderable;

void SingleResolveRenderer::Init(
    std::function<std::unique_ptr<Scene>()> scene_factory,
    const Renderable::RenderingContext& context) {
  scene_ = scene_factory();
  context_ = context;
}

void SingleResolveRenderer::SetComponent(
    std::unique_ptr<const Component> root) {
  CHECK_NOTNULL(scene_);
  scene_->Clear();
  root_component_ = std::move(root);

  // Get the node set and create a root node with the opaque node and the
  // transparent nodes as children.
  Renderable::NodeSet node_set =
      root_component_->GetRenderable()->GetRenderNodes(context_);
  ion::gfx::NodePtr root_node(new ion::gfx::Node);
  root_node->AddChild(node_set.opaque);
  root_node->AddChild(node_set.transparent);
  scene_->AddNode(root_node);
}

void SingleResolveRenderer::RenderFrame(
    const Renderable::ViewState& view_state) {
  root_component_->GetRenderable()->Update(view_state);
  context_.renderer->DrawScene(scene_->GetRoot());
}

}  // namespace viewer
}  // namespace seurat
