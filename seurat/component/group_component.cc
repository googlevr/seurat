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

#include "seurat/component/group_component.h"

#include "ion/gfx/node.h"
#include "seurat/base/structured_io.h"
#include "seurat/component/component.h"
#include "seurat/component/renderable.h"

namespace seurat {
namespace component {

namespace {

class GroupRenderable : public Renderable {
 public:
  explicit GroupRenderable(const GroupComponent* group_component)
      : group_component_(group_component) {
    CHECK_NOTNULL(group_component);
  }

  ~GroupRenderable() override = default;

  Properties GetProperties() const override {
    Properties properties;

    for (const auto& child : group_component_->GetChildren()) {
      Properties child_props = child->GetRenderable()->GetProperties();

      // Merge child_props into properties.
      properties.bounding_box.ExtendByRange(child_props.bounding_box);
    }

    return properties;
  }

  NodeSet GetRenderNodes(const RenderingContext& context) override {
    NodeSet node_set;
    std::vector<ion::gfx::NodePtr*> all_nodes = node_set.AllNodes();
    for (const auto& child : group_component_->GetChildren()) {
      NodeSet child_nodes = child->GetRenderable()->GetRenderNodes(context);
      std::vector<ion::gfx::NodePtr*> all_child_nodes = child_nodes.AllNodes();

      // Merge all_child_nodes into all_nodes
      for (int i = 0; i < all_nodes.size(); ++i) {
        if (!all_child_nodes[i]->Get()) {
          // If the child's i'th node is null, skip it.
          continue;
        }

        // If this is the first time a child has had a nonnull i'th node,
        // allocate a new node to store it.
        if (!all_nodes[i]->Get()) {
          *all_nodes[i] = ion::gfx::NodePtr(new ion::gfx::Node);
          (*all_nodes[i])->SetLabel(group_component_->GetLabel());
        }
        all_nodes[i]->Get()->AddChild(*all_child_nodes[i]);
      }
    }
    return node_set;
  }

  void Update(const ViewState& view_state) override {
    for (const auto& child : group_component_->GetChildren()) {
      child->GetRenderable()->Update(view_state);
    }
  }

 private:
  const GroupComponent* group_component_;
};

}  // namespace

REGISTER_SEURAT_COMPONENT(GroupComponent);

GroupComponent::GroupComponent(
    std::string label, std::vector<std::unique_ptr<const Component>> children)
    : Component(std::move(label)),
      children_(std::move(children)),
      renderable_(new GroupRenderable(this)) {}

// static
std::unique_ptr<const Component> GroupComponent::Create(
    std::string label, base::StructureSource* source) {
  std::vector<std::unique_ptr<const Component>> children;
  uint64 child_count = source->ReadPod<uint64>();
  children.reserve(child_count);
  for (int i = 0; i < child_count; ++i) {
    children.emplace_back(Component::Create(source));
  }
  return std::unique_ptr<const Component>(
      new GroupComponent(std::move(label), std::move(children)));
}

void GroupComponent::WriteInternal(base::StructureSink* sink) const {
  sink->WritePod<uint64>(children_.size());
  for (const auto& child : children_) {
    child->Write(sink);
  }
}

bool GroupComponent::operator==(const Component& other) const {
  if (!Component::IsComponentEqualTo(other)) {
    return false;
  }
  const GroupComponent& ogroup = static_cast<const GroupComponent&>(other);
  if (children_.size() != ogroup.children_.size()) {
    return false;
  }
  for (int i = 0; i < children_.size(); ++i) {
    if (*children_[i] != *ogroup.children_[i]) {
      return false;
    }
  }
  return true;
}

// static
std::unique_ptr<const GroupComponent> GroupComponent::CreateFromComponents(
    std::string label, std::unique_ptr<const Component> first,
    std::unique_ptr<const Component> second) {
  std::vector<std::unique_ptr<const Component>> children(2);
  children[0] = std::move(first);
  children[1] = std::move(second);
  std::unique_ptr<const GroupComponent> group_component(
      new component::GroupComponent(std::move(label), std::move(children)));
  return group_component;
}

}  // namespace component
}  // namespace seurat
