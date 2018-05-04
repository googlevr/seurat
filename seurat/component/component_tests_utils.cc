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

#include "seurat/component/component_tests_utils.h"

#include "gtest/gtest.h"
#include "absl/strings/string_view.h"
#include "seurat/base/structured_io.h"

namespace seurat {
namespace component {

using base::StructureSink;
using base::StructureSource;
using ion::math::Matrix4f;

REGISTER_SEURAT_COMPONENT(FakeComponent);

FakeComponent::FakeComponent(std::string label,
                             std::unique_ptr<Renderable> renderable)
    : Component(std::move(label)), renderable_(std::move(renderable)) {}

// static
std::unique_ptr<const Component> FakeComponent::Create(
    std::string label, StructureSource* source) {
  EXPECT_EQ(3, source->ReadPod<uint8>());
  EXPECT_EQ(1, source->ReadPod<uint8>());
  EXPECT_EQ(4, source->ReadPod<uint8>());
  EXPECT_EQ(1, source->ReadPod<uint8>());
  EXPECT_EQ(5, source->ReadPod<uint8>());
  return std::unique_ptr<const Component>(
      new FakeComponent(std::move(label), nullptr));
}

void FakeComponent::WriteInternal(StructureSink* sink) const {
  sink->WritePod<uint8>(3);
  sink->WritePod<uint8>(1);
  sink->WritePod<uint8>(4);
  sink->WritePod<uint8>(1);
  sink->WritePod<uint8>(5);
}

bool FakeComponent::operator==(const Component& other) const {
  return IsComponentEqualTo(other);
}

Renderable* FakeComponent::GetRenderable() const { return renderable_.get(); }

FakeRenderable::FakeRenderable()
    : nodes_(),
      props_(),
      num_updates_(0),
      left_clip_from_component_matrix_(Matrix4f::Identity()),
      right_clip_from_component_matrix_(Matrix4f::Identity()) {}

FakeRenderable::FakeRenderable(NodeSet nodes)
    : nodes_(nodes),
      props_(),
      num_updates_(0),
      left_clip_from_component_matrix_(Matrix4f::Identity()),
      right_clip_from_component_matrix_(Matrix4f::Identity()) {}

FakeRenderable::FakeRenderable(Properties props)
    : nodes_(),
      props_(props),
      num_updates_(0),
      left_clip_from_component_matrix_(Matrix4f::Identity()),
      right_clip_from_component_matrix_(Matrix4f::Identity()) {}

Renderable::Properties FakeRenderable::GetProperties() const { return props_; }

Renderable::NodeSet FakeRenderable::GetRenderNodes(
    const RenderingContext& context) {
  return nodes_;
}

void FakeRenderable::Update(const ViewState& view_state) {
  num_updates_++;
  left_clip_from_component_matrix_ = view_state.left_clip_from_component_matrix;
  right_clip_from_component_matrix_ =
      view_state.right_clip_from_component_matrix;
}

int FakeRenderable::GetNumUpdates() const { return num_updates_; }

const Matrix4f& FakeRenderable::GetLeftClipFromComponentMatrix() const {
  return left_clip_from_component_matrix_;
}

const Matrix4f& FakeRenderable::GetRightClipFromComponentMatrix() const {
  return right_clip_from_component_matrix_;
}

// static
void ComponentCommonTests::TestIo(const Component& a, const Component& b) {
  ASSERT_NE(a, b);

  std::unique_ptr<Renderable> fake_renderable(new FakeRenderable());
  std::unique_ptr<FakeComponent> fake_component(
      new FakeComponent("fake", std::move(fake_renderable)));
  EXPECT_EQ(a, a);
  EXPECT_EQ(b, b);
  EXPECT_NE(*fake_component, a);
  EXPECT_NE(a, *fake_component);
  EXPECT_NE(b, *fake_component);

  std::string str;
  base::StringByteSink byte_sink(&str);
  StructureSink component_sink(&byte_sink);

  a.Write(&component_sink);
  b.Write(&component_sink);
  fake_component->Write(&component_sink);

  size_t num_bytes_written = str.size();

  base::ArrayByteSource byte_source(
      absl::string_view(str.data(), num_bytes_written));
  StructureSource component_source(&byte_source);

  const std::unique_ptr<const Component> component_a_read =
      Component::Create(&component_source);
  const std::unique_ptr<const Component> component_b_read =
      Component::Create(&component_source);
  const std::unique_ptr<const Component> fake_component_read =
      Component::Create(&component_source);

  EXPECT_EQ(a, *component_a_read);
  EXPECT_EQ(*component_a_read, a);
  EXPECT_EQ(b, *component_b_read);
  EXPECT_EQ(*component_b_read, b);
  EXPECT_EQ(*fake_component, *fake_component_read);

  EXPECT_NE(a, *component_b_read);
  EXPECT_NE(*component_b_read, a);
  EXPECT_NE(b, *component_a_read);
  EXPECT_NE(*component_a_read, b);
  EXPECT_NE(*fake_component, *component_a_read);
  EXPECT_NE(*component_a_read, *fake_component);
  EXPECT_NE(*fake_component, *component_b_read);
  EXPECT_NE(*component_b_read, *fake_component);
}

}  // namespace component
}  // namespace seurat
