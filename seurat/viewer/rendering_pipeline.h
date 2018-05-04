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

#ifndef VR_SEURAT_VIEWER_RENDERING_PIPELINE_H_
#define VR_SEURAT_VIEWER_RENDERING_PIPELINE_H_

#include <functional>
#include <memory>

#include "seurat/component/component.h"
#include "seurat/component/renderable.h"
#include "seurat/viewer/scene.h"

namespace seurat {
namespace viewer {

class Renderer {
 public:
  virtual ~Renderer() = default;

  // Initializes the renderer with the given |context|. This is usually invoked
  // once, during the initialization stage of the application, or when the
  // context changes (e.g. window resize).
  virtual void Init(std::function<std::unique_ptr<Scene>()> scene_factory,
                    const component::Renderable::RenderingContext& context) = 0;

  // Returns a pointer to the root component.
  virtual const component::Component* GetComponent() const = 0;

  // Initializes or replaces (in the case of keyframe animation) the current
  // root Component to render.
  virtual void SetComponent(
      std::unique_ptr<const component::Component> root) = 0;

  // Renders a single frame.
  virtual void RenderFrame(
      const component::Renderable::ViewState& view_state) = 0;

  // Returns the root of the render nodes.
  virtual const ion::gfx::NodePtr& GetSceneRoot() = 0;
};

// A renderer which renders on only a single framebuffer resolve (excluding any
// VR lens-distortion pass).
//
// Note that all Nodes in the seurat::Renderable::NodeSet are *ignored*
// except:
//  * opaque
//  * transparent
class SingleResolveRenderer : public Renderer {
 public:
  SingleResolveRenderer() = default;
  ~SingleResolveRenderer() override = default;

  void Init(std::function<std::unique_ptr<Scene>()> scene_factory,
            const component::Renderable::RenderingContext& context) override;

  const component::Component* GetComponent() const override {
    return root_component_.get();
  }

  void SetComponent(std::unique_ptr<const component::Component> root) override;

  void RenderFrame(const component::Renderable::ViewState& view_state) override;

  const ion::gfx::NodePtr& GetSceneRoot() override { return scene_->GetRoot(); }

 private:
  // The context for this renderer. Only valid after Init was called.
  component::Renderable::RenderingContext context_;

  // The scene that holds the root component's render nodes.
  std::unique_ptr<Scene> scene_;

  // The root node of the component hierarchy.
  std::unique_ptr<const component::Component> root_component_;
};

}  // namespace viewer
}  // namespace seurat

#endif  // VR_SEURAT_VIEWER_RENDERING_PIPELINE_H_
