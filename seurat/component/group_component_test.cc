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

#include <memory>

#include "gtest/gtest.h"
#include "seurat/component/component.h"
#include "seurat/component/component_tests_utils.h"
#include "seurat/component/renderable.h"

namespace seurat {
namespace component {
namespace {

// Returns a pre-order traversal of the given ion node tree into |all_nodes|.
void PreorderTraversal(ion::gfx::NodePtr root,
                       std::vector<ion::gfx::NodePtr>* all_nodes) {
  if (!root.Get()) return;
  all_nodes->push_back(root);
  for (const auto& child : root->GetChildren()) {
    PreorderTraversal(child, all_nodes);
  }
}

TEST(GroupComponentTest, TestIo) {
  std::unique_ptr<Component> child1(new GroupComponent());
  std::unique_ptr<Component> child2(new GroupComponent());
  std::unique_ptr<Component> child3(new GroupComponent());
  std::vector<std::unique_ptr<const Component>> components_a;
  components_a.push_back(std::move(child1));

  std::vector<std::unique_ptr<const Component>> components_b;
  components_b.push_back(std::move(child2));
  components_b.push_back(std::move(child3));

  GroupComponent a(std::move(components_a));
  GroupComponent b(std::move(components_b));

  ComponentCommonTests::TestIo(a, b);
}

TEST(GroupComponentTest, RenderablePropertiesShouldBeCombined) {
  auto props1 = Renderable::Properties{
      ion::math::Range3f({2.0f, 4.0f, 6.0f}, {8.0f, 10.0f, 12.0f})};
  auto props2 = Renderable::Properties{
      ion::math::Range3f({4.0f, 2.0f, 4.0f}, {6.0f, 8.0f, 14.0f})};

  auto renderable1 =
      std::unique_ptr<FakeRenderable>(new FakeRenderable(props1));
  auto renderable2 =
      std::unique_ptr<FakeRenderable>(new FakeRenderable(props2));

  std::vector<std::unique_ptr<const Component>> children;
  children.push_back(std::unique_ptr<Component>(
      new FakeComponent("1", std::move(renderable1))));
  children.push_back(std::unique_ptr<Component>(
      new FakeComponent("2", std::move(renderable2))));
  GroupComponent group(std::move(children));

  auto combined_props = group.GetRenderable()->GetProperties();

  EXPECT_EQ(ion::math::Range3f({2.0f, 2.0f, 4.0f}, {8.0f, 10.0f, 14.0f}),
            combined_props.bounding_box);
}

TEST(GroupComponentTest, ChildRenderablesShouldBeUpdated) {
  auto renderable1 = std::unique_ptr<FakeRenderable>(new FakeRenderable());
  auto renderable2 = std::unique_ptr<FakeRenderable>(new FakeRenderable());

  FakeRenderable* renderable1p = renderable1.get();
  FakeRenderable* renderable2p = renderable2.get();

  std::vector<std::unique_ptr<const Component>> children;
  children.push_back(std::unique_ptr<Component>(
      new FakeComponent("1", std::move(renderable1))));
  children.push_back(std::unique_ptr<Component>(
      new FakeComponent("2", std::move(renderable2))));
  GroupComponent group(std::move(children));

  auto* group_renderable = group.GetRenderable();

  group_renderable->Update({});
  EXPECT_EQ(1, renderable1p->GetNumUpdates());
  EXPECT_EQ(1, renderable2p->GetNumUpdates());

  group_renderable->Update({});
  EXPECT_EQ(2, renderable1p->GetNumUpdates());
  EXPECT_EQ(2, renderable2p->GetNumUpdates());
}

TEST(GroupComponentTest, CreateFromComponents) {
  auto renderable1 = std::unique_ptr<FakeRenderable>(new FakeRenderable());
  auto renderable2 = std::unique_ptr<FakeRenderable>(new FakeRenderable());

  std::unique_ptr<Component> first_child(
      new FakeComponent("1", std::move(renderable1)));
  std::unique_ptr<Component> second_child(
      new FakeComponent("2", std::move(renderable2)));
  const Component* first_child_pointer = first_child.get();
  const Component* second_child_pointer = second_child.get();

  std::unique_ptr<const GroupComponent> group =
      GroupComponent::CreateFromComponents("group", std::move(first_child),
                                           std::move(second_child));
  EXPECT_EQ(group->GetChildren().size(), 2);
  EXPECT_EQ(group->GetChildren()[0].get(), first_child_pointer);
  EXPECT_EQ(group->GetChildren()[1].get(), second_child_pointer);
}

TEST(GroupComponentTest, EmptyGroupShouldHaveEmptyNodeSet) {
  std::vector<std::unique_ptr<const Component>> children;
  GroupComponent empty_group(std::move(children));

  Renderable::NodeSet render_nodes =
      empty_group.GetRenderable()->GetRenderNodes({});
  EXPECT_EQ(nullptr, render_nodes.opaque.Get());
  EXPECT_EQ(nullptr, render_nodes.transparent.Get());
}

TEST(GroupComponentTest, ChildIonNodesShouldBeMergedInOrder) {
  // Construct a GroupComponent with the following children:
  //  * child1 - Renders a opaque pass
  //  * child2 - Renders nothing
  //  * child3 - Renders an opaque pass and a transparent pass

  Renderable::NodeSet nodes_1;
  Renderable::NodeSet nodes_2;
  Renderable::NodeSet nodes_3;

  nodes_1.opaque = ion::gfx::NodePtr(new ion::gfx::Node);
  nodes_1.opaque->SetLabel("Node1 Opaque");
  nodes_3.opaque = ion::gfx::NodePtr(new ion::gfx::Node);
  nodes_3.opaque->SetLabel("Node3 Opaque");
  nodes_3.transparent = ion::gfx::NodePtr(new ion::gfx::Node);
  nodes_3.transparent->SetLabel("Node3 Transparent");

  std::vector<std::unique_ptr<const Component>> children;
  {
    std::unique_ptr<Component> child1(new FakeComponent(
        "1", std::unique_ptr<Renderable>(new FakeRenderable(nodes_1))));
    std::unique_ptr<Component> child2(new FakeComponent(
        "2", std::unique_ptr<Renderable>(new FakeRenderable(nodes_2))));
    std::unique_ptr<Component> child3(new FakeComponent(
        "3", std::unique_ptr<Renderable>(new FakeRenderable(nodes_3))));

    children.push_back(std::move(child1));
    children.push_back(std::move(child2));
    children.push_back(std::move(child3));
  }
  GroupComponent group(std::move(children));

  Renderable::NodeSet group_nodes = group.GetRenderable()->GetRenderNodes({});

  EXPECT_NE(nullptr, group_nodes.opaque.Get());
  EXPECT_NE(nullptr, group_nodes.transparent.Get());

  std::vector<ion::gfx::NodePtr> all_opaque_nodes;
  PreorderTraversal(group_nodes.opaque, &all_opaque_nodes);

  std::vector<ion::gfx::NodePtr> all_transparent_nodes;
  PreorderTraversal(group_nodes.transparent, &all_transparent_nodes);

  // Opaque nodes must be present.
  EXPECT_EQ(1, std::count(all_opaque_nodes.begin(), all_opaque_nodes.end(),
                          nodes_1.opaque));
  EXPECT_EQ(1, std::count(all_opaque_nodes.begin(), all_opaque_nodes.end(),
                          nodes_3.opaque));

  // Opaque should not be found in the color pass.
  EXPECT_EQ(0, std::count(all_transparent_nodes.begin(),
                          all_transparent_nodes.end(), nodes_1.opaque));
  EXPECT_EQ(0, std::count(all_transparent_nodes.begin(),
                          all_transparent_nodes.end(), nodes_3.opaque));

  // Opaque nodes must be in order.
  int child_1_opaque_index = std::find(all_opaque_nodes.begin(),
                                       all_opaque_nodes.end(), nodes_1.opaque) -
                             all_opaque_nodes.begin();
  int child_3_opaque_index = std::find(all_opaque_nodes.begin(),
                                       all_opaque_nodes.end(), nodes_3.opaque) -
                             all_opaque_nodes.begin();
  EXPECT_GT(child_3_opaque_index, child_1_opaque_index);

  EXPECT_EQ(1, std::count(all_transparent_nodes.begin(),
                          all_transparent_nodes.end(), nodes_3.transparent));
}

}  // namespace
}  // namespace component
}  // namespace seurat
