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

#ifndef VR_SEURAT_COMPONENT_GROUP_COMPONENT_H_
#define VR_SEURAT_COMPONENT_GROUP_COMPONENT_H_

#include <memory>
#include <vector>

#include "seurat/base/structured_io.h"
#include "seurat/component/component.h"
#include "seurat/component/renderable.h"

namespace seurat {
namespace component {

// Joins a list of components which are all rendered together, in order.
class GroupComponent : public Component {
 public:
  GroupComponent(std::string label,
                 std::vector<std::unique_ptr<const Component>> children);
  explicit GroupComponent(std::string label)
      : GroupComponent(std::move(label), {}) {}
  GroupComponent() : GroupComponent({}, {}) {}
  explicit GroupComponent(
      std::vector<std::unique_ptr<const Component>> children)
      : GroupComponent({}, std::move(children)) {}
  ~GroupComponent() override = default;

  static std::unique_ptr<const Component> Create(std::string label,
                                                 base::StructureSource* source);
  bool operator==(const Component& other) const override;
  Renderable* GetRenderable() const override { return renderable_.get(); }

  const std::vector<std::unique_ptr<const Component>>& GetChildren() const {
    return children_;
  }

  // Combines two input components |first| and |second| in order with a
  // GroupComponent, labelling the combination with |label|, and returns the
  // composite component.
  static std::unique_ptr<const GroupComponent> CreateFromComponents(
      std::string label, std::unique_ptr<const Component> first,
      std::unique_ptr<const Component> second);

 private:
  void WriteInternal(base::StructureSink* sink) const override;

  std::vector<std::unique_ptr<const Component>> children_;
  std::unique_ptr<Renderable> renderable_;

  DECLARE_SEURAT_COMPONENT(GroupComponent);
};

}  // namespace component
}  // namespace seurat

#endif  // VR_SEURAT_COMPONENT_GROUP_COMPONENT_H_
