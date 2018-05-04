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

#ifndef VR_SEURAT_COMPONENT_TRANSFORM_COMPONENT_H_
#define VR_SEURAT_COMPONENT_TRANSFORM_COMPONENT_H_

#include <memory>
#include <vector>

#include "seurat/component/component.h"
#include "seurat/component/renderable.h"

namespace seurat {
namespace component {

// Transforms its children.
class TransformComponent : public Component {
 public:
  explicit TransformComponent(std::string label);
  TransformComponent() : TransformComponent("") {}
  TransformComponent(std::string label,
                     const ion::math::Matrix4f& parent_from_child_matrix,
                     std::unique_ptr<const Component> child);
  ~TransformComponent() override = default;

  static std::unique_ptr<const Component> Create(std::string label,
                                                 base::StructureSource* source);
  bool operator==(const Component& other) const override;
  Renderable* GetRenderable() const override { return renderable_.get(); }

  const Component* GetChild() const { return child_.get(); }
  const ion::math::Matrix4f& GetParentFromChildMatrix() const {
    return parent_from_child_matrix_;
  }

  // Wraps an input |component| with a transformation |matrix| via
  // TransformComponent, labelling the transformation node with
  // |label|, and returns the composite transformed component.
  static std::unique_ptr<const TransformComponent> CreateFromComponent(
      std::string label, const ion::math::Matrix4f& matrix,
      std::unique_ptr<const Component> component);

 private:
  void WriteInternal(base::StructureSink* sink) const override;

  // Transformation from the child's coordinate space to the parent's coordinate
  // space.
  ion::math::Matrix4f parent_from_child_matrix_;
  std::unique_ptr<const Component> child_;
  std::unique_ptr<Renderable> renderable_;

  DECLARE_SEURAT_COMPONENT(TransformComponent);
};

}  // namespace component
}  // namespace seurat

#endif  // VR_SEURAT_COMPONENT_TRANSFORM_COMPONENT_H_
