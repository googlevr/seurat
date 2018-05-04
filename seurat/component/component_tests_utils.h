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

#ifndef VR_SEURAT_COMPONENT_COMPONENT_TESTS_UTILS_H_
#define VR_SEURAT_COMPONENT_COMPONENT_TESTS_UTILS_H_

#include <array>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <type_traits>

#include "ion/math/matrix.h"
#include "ion/math/vector.h"
#include "seurat/base/util.h"
#include "seurat/component/component.h"
#include "seurat/component/renderable.h"

namespace seurat {
namespace component {

// A Component for testing which can be initialized with any renderable.
class FakeComponent : public Component {
 public:
  FakeComponent(std::string label, std::unique_ptr<Renderable> renderable);
  explicit FakeComponent(std::unique_ptr<Renderable> renderable)
      : FakeComponent({}, std::move(renderable)) {}
  ~FakeComponent() override = default;

  // Component implementation.
  static std::unique_ptr<const Component> Create(std::string label,
                                                 base::StructureSource* source);
  void WriteInternal(base::StructureSink* sink) const override;
  bool operator==(const Component& other) const override;
  Renderable* GetRenderable() const override;

 private:
  const std::unique_ptr<Renderable> renderable_;

  DECLARE_SEURAT_COMPONENT(FakeComponent);
};

// A Renderable for testing which can be initialized with any NodeSet and
// Properties.
class FakeRenderable : public Renderable {
 public:
  FakeRenderable();
  explicit FakeRenderable(NodeSet nodes);
  explicit FakeRenderable(Properties props);
  ~FakeRenderable() override = default;

  // Renderable implementation.
  Properties GetProperties() const override;
  NodeSet GetRenderNodes(const RenderingContext& context) override;
  void Update(const ViewState& view_state) override;

  int GetNumUpdates() const;
  const ion::math::Matrix4f& GetLeftClipFromComponentMatrix() const;
  const ion::math::Matrix4f& GetRightClipFromComponentMatrix() const;

 private:
  const NodeSet nodes_;
  const Properties props_;
  int num_updates_;
  // The 'cached' clip-from-component matrices. These are set by Update() to the
  // matrices of the incoming view_state.
  ion::math::Matrix4f left_clip_from_component_matrix_;
  ion::math::Matrix4f right_clip_from_component_matrix_;
};

// Contains tests which should apply to all component implementations.
class ComponentCommonTests {
 public:
  // Tests Read & Write methods for two components.  The given example
  // instances, |a| and |b|, must be different.
  static void TestIo(const Component& a, const Component& b);
};

}  // namespace component
}  // namespace seurat

#endif  // VR_SEURAT_COMPONENT_COMPONENT_TESTS_UTILS_H_
