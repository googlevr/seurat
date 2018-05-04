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

#include "seurat/component/transform_component.h"

#include "ion/gfx/node.h"
#include "ion/math/matrix.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/base/structured_io.h"
#include "seurat/component/component.h"
#include "seurat/component/renderable.h"

using ion::math::Matrix4f;

namespace seurat {
namespace component {

namespace {

class TransformRenderable : public Renderable {
 public:
  explicit TransformRenderable(const TransformComponent* transform_component)
      : transform_component_(transform_component) {
    CHECK_NOTNULL(transform_component);
  }

  ~TransformRenderable() override = default;

  Properties GetProperties() const override {
    Properties properties;
    const Component* child = transform_component_->GetChild();
    if (child) {
      // Transform the child's AABB.
      const Matrix4f parent_from_child =
          transform_component_->GetParentFromChildMatrix();
      properties.bounding_box = base::ProjectAABB(
          parent_from_child,
          child->GetRenderable()->GetProperties().bounding_box);
    }
    return properties;
  }

  NodeSet GetRenderNodes(const RenderingContext& context) override {
    const Component* child = transform_component_->GetChild();
    if (child) {
      return child->GetRenderable()->GetRenderNodes(context);
    } else {
      return NodeSet();
    }
  }

  void Update(const ViewState& view_state) override {
    const Component* child = transform_component_->GetChild();
    if (!child) {
      return;
    }
    ViewState transformed_view_state = view_state;
    const Matrix4f parent_from_child_matrix =
        transform_component_->GetParentFromChildMatrix();
    transformed_view_state.left_clip_from_component_matrix =
        transformed_view_state.left_clip_from_component_matrix *
        parent_from_child_matrix;
    transformed_view_state.right_clip_from_component_matrix =
        transformed_view_state.right_clip_from_component_matrix *
        parent_from_child_matrix;
    transformed_view_state.left_eye_from_component_matrix =
        transformed_view_state.left_eye_from_component_matrix *
        parent_from_child_matrix;
    transformed_view_state.right_eye_from_component_matrix =
        transformed_view_state.right_eye_from_component_matrix *
        parent_from_child_matrix;
    child->GetRenderable()->Update(transformed_view_state);
  }

 private:
  const TransformComponent* transform_component_;
};

}  // namespace

REGISTER_SEURAT_COMPONENT(TransformComponent);

TransformComponent::TransformComponent(std::string label)
    : Component(std::move(label)),
      parent_from_child_matrix_(Matrix4f::Identity()),
      child_(nullptr),
      renderable_(new TransformRenderable(this)) {}

TransformComponent::TransformComponent(std::string label,
                                       const Matrix4f& parent_from_child_matrix,
                                       std::unique_ptr<const Component> child)
    : Component(std::move(label)),
      parent_from_child_matrix_(parent_from_child_matrix),
      child_(std::move(child)),
      renderable_(new TransformRenderable(this)) {}

// static
std::unique_ptr<const Component> TransformComponent::Create(
    std::string label, base::StructureSource* source) {
  const Matrix4f parent_from_child_matrix = source->ReadMatrix<Matrix4f>();
  std::unique_ptr<const Component> child;
  const bool has_child = source->ReadPod<bool>();
  if (has_child) {
    child = Component::Create(source);
  }
  return std::unique_ptr<const Component>(new TransformComponent(
      std::move(label), parent_from_child_matrix, std::move(child)));
}

void TransformComponent::WriteInternal(base::StructureSink* sink) const {
  sink->WriteMatrix(parent_from_child_matrix_);
  sink->WritePod<bool>(child_ ? true : false);
  if (child_) {
    child_->Write(sink);
  }
}

bool TransformComponent::operator==(const Component& other) const {
  if (!IsComponentEqualTo(other)) {
    return false;
  }
  const TransformComponent& other_transform =
      static_cast<const TransformComponent&>(other);
  if (other_transform.parent_from_child_matrix_ != parent_from_child_matrix_) {
    return false;
  }
  if (child_ && (*child_ != *other_transform.child_)) {
    return false;
  }
  return true;
}

// static
std::unique_ptr<const TransformComponent>
TransformComponent::CreateFromComponent(
    std::string transformed_label, const ion::math::Matrix4f& matrix,
    std::unique_ptr<const Component> component) {
  std::unique_ptr<const TransformComponent> transformed_component(
      new TransformComponent(std::move(transformed_label), matrix,
                             std::move(component)));
  return transformed_component;
}

}  // namespace component
}  // namespace seurat
