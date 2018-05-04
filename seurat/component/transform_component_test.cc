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

#include <memory>

#include "ion/math/matrix.h"
#include "ion/math/transformutils.h"
#include "gtest/gtest.h"
#include "seurat/component/component_tests_utils.h"
#include "seurat/component/group_component.h"

using ion::math::Matrix4f;
using ion::math::Vector3f;

namespace seurat {
namespace component {
namespace {

TEST(TransformComponentTest, TestIo) {
  Matrix4f matrix_a = ion::math::TranslationMatrix(Vector3f(1.0f, -2.0f, 3.0f));
  std::unique_ptr<Component> child_a(new GroupComponent());
  TransformComponent a("a", matrix_a, std::move(child_a));

  Matrix4f matrix_b = ion::math::ScaleMatrixH(Vector3f(-0.5f, 1.5f, 2.0f));
  std::unique_ptr<Component> child_b(new GroupComponent());
  TransformComponent b("b", matrix_b, std::move(child_b));

  ComponentCommonTests::TestIo(a, b);
}

TEST(TransformComponentTest, RenderablePropertiesArePassedThrough) {
  auto props = Renderable::Properties{
      ion::math::Range3f({2.0f, 4.0f, 6.0f}, {8.0f, 10.0f, 12.0f})};
  auto renderable = std::unique_ptr<FakeRenderable>(new FakeRenderable(props));
  Matrix4f matrix = ion::math::TranslationMatrix(Vector3f(1.0f, -2.0f, 3.0f));
  std::unique_ptr<const Component> child(
      new FakeComponent("child", std::move(renderable)));
  TransformComponent transform("transformer", matrix, std::move(child));

  auto actual_props = transform.GetRenderable()->GetProperties();

  EXPECT_EQ(ion::math::Range3f({3.0f, 2.0f, 9.0f}, {9.0f, 8.0f, 15.0f}),
            actual_props.bounding_box);
}

TEST(TransformComponentTest, TransformWithNoChildShouldHaveEmptyNodeSet) {
  TransformComponent leaf_transform;
  Renderable::NodeSet render_nodes =
      leaf_transform.GetRenderable()->GetRenderNodes({});
  EXPECT_EQ(nullptr, render_nodes.opaque.Get());
  EXPECT_EQ(nullptr, render_nodes.transparent.Get());
}

TEST(TransformComponentTest, NodeSetIsPassedThrough) {
  Renderable::NodeSet nodes;
  nodes.opaque = ion::gfx::NodePtr(new ion::gfx::Node);
  nodes.opaque->SetLabel("opaque");
  nodes.transparent = ion::gfx::NodePtr(new ion::gfx::Node);
  nodes.transparent->SetLabel("transparent");
  std::unique_ptr<Component> child(new FakeComponent(
      "child", std::unique_ptr<Renderable>(new FakeRenderable(nodes))));

  Matrix4f matrix = Matrix4f::Identity();
  TransformComponent transform("identity", matrix, std::move(child));
  Renderable::NodeSet render_nodes =
      transform.GetRenderable()->GetRenderNodes({});
  EXPECT_EQ(nodes.opaque.Get(), render_nodes.opaque.Get());
  EXPECT_EQ(nodes.transparent.Get(), render_nodes.transparent.Get());
}

TEST(TransformComponentTest, CreateFromComponent) {
  auto props = Renderable::Properties{
      ion::math::Range3f({0.0f, 0.0f, 0.0f}, {1.0f, 1.0f, 1.0f})};
  auto renderable = std::unique_ptr<FakeRenderable>(new FakeRenderable(props));
  Matrix4f matrix = ion::math::TranslationMatrix(Vector3f(1.0f, -2.0f, 3.0f));
  std::unique_ptr<const Component> child(
      new FakeComponent("child", std::move(renderable)));
  const Component* child_pointer = child.get();
  std::unique_ptr<const TransformComponent> transform =
      TransformComponent::CreateFromComponent("translate", matrix,
                                              std::move(child));
  EXPECT_EQ(transform->GetChild(), child_pointer);
}

TEST(TransformComponentTest, MatricesArePropagatedDown) {
  auto props = Renderable::Properties{
      ion::math::Range3f({0.0f, 0.0f, 0.0f}, {1.0f, 1.0f, 1.0f})};
  auto renderable = std::unique_ptr<FakeRenderable>(new FakeRenderable(props));
  FakeRenderable* fake_renderable = renderable.get();
  Matrix4f matrix = ion::math::TranslationMatrix(Vector3f(1.0f, -2.0f, 3.0f));
  std::unique_ptr<const Component> child(
      new FakeComponent("child", std::move(renderable)));
  TransformComponent transform("translate", matrix, std::move(child));
  EXPECT_EQ(0, fake_renderable->GetNumUpdates());
  EXPECT_EQ(Matrix4f::Identity(),
            fake_renderable->GetLeftClipFromComponentMatrix());
  EXPECT_EQ(Matrix4f::Identity(),
            fake_renderable->GetRightClipFromComponentMatrix());

  // Call Update() and expect the matrix to be propagated down the hierarchy.
  Renderable::ViewState view_state;
  view_state.left_clip_from_component_matrix = Matrix4f::Identity();
  view_state.right_clip_from_component_matrix = -1.0f * Matrix4f::Identity();
  transform.GetRenderable()->Update(view_state);
  EXPECT_EQ(1, fake_renderable->GetNumUpdates());
  EXPECT_EQ(matrix, fake_renderable->GetLeftClipFromComponentMatrix());
  EXPECT_EQ(-1.0f * matrix, fake_renderable->GetRightClipFromComponentMatrix());
}

}  // namespace
}  // namespace component
}  // namespace seurat
