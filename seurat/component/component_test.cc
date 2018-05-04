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

#include "seurat/component/component.h"

#include <memory>

#include "ion/math/matrix.h"
#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "absl/strings/string_view.h"
#include "seurat/base/bytestream.h"
#include "seurat/base/structured_io.h"
#include "seurat/component/component_tests_utils.h"
#include "seurat/component/renderable.h"

namespace seurat {
namespace component {
namespace {

using base::StructureSink;
using base::StructureSource;
using ion::math::Matrix4f;
using ion::math::Point2i;
using ion::math::Vector3i;

class ComponentA : public Component {
 public:
  explicit ComponentA(std::string label = {})
      : Component(std::move(label)), renderable_(new FakeRenderable()) {}
  ~ComponentA() override = default;

  int GetIValue() const { return ivalue_; }
  float GetFValue() const { return fvalue_; }
  void SetIValue(int value) { ivalue_ = value; }
  void SetFValue(float value) { fvalue_ = value; }

  // Component implementation.
  static std::unique_ptr<const Component> Create(std::string label,
                                                 StructureSource* source) {
    std::unique_ptr<ComponentA> component(new ComponentA(std::move(label)));
    component->SetIValue(source->ReadPod<int>());
    component->SetFValue(source->ReadPod<float>());
    return std::move(component);
  }

  void WriteInternal(StructureSink* sink) const override {
    sink->WritePod<int>(ivalue_);
    sink->WritePod<float>(fvalue_);
  }

  bool operator==(const Component& other) const override {
    if (!IsComponentEqualTo(other)) {
      return false;
    }
    const ComponentA& other_a = static_cast<const ComponentA&>(other);
    return ivalue_ == other_a.ivalue_ && fvalue_ == other_a.fvalue_;
  }

  Renderable* GetRenderable() const override { return renderable_.get(); }

 private:
  const std::unique_ptr<FakeRenderable> renderable_;
  int ivalue_;
  float fvalue_;

  DECLARE_SEURAT_COMPONENT(ComponentA);
};
REGISTER_SEURAT_COMPONENT(ComponentA);

class ComponentB : public Component {
 public:
  explicit ComponentB(std::string label = {})
      : Component(std::move(label)),
        renderable_(new FakeRenderable),
        vec_({1, 2, 3}),
        point_({4, 5}),
        mat_(Matrix4f::Identity()) {}
  ~ComponentB() override = default;

  Vector3i GetVec() const { return vec_; }
  Point2i GetPoint() const { return point_; }
  Matrix4f GetMatrix() const { return mat_; }
  const Component* GetChild() const { return child_.get(); }
  void SetChild(std::unique_ptr<const Component> child) {
    child_ = std::move(child);
  }

  // Component implementation.
  static std::unique_ptr<const Component> Create(std::string label,
                                                 StructureSource* source) {
    std::unique_ptr<ComponentB> component(new ComponentB(std::move(label)));
    component->vec_ = source->ReadVector<Vector3i>();
    component->point_ = source->ReadPoint<Point2i>();
    component->mat_ = source->ReadMatrix<Matrix4f>();
    const bool has_child = source->ReadPod<bool>();
    if (has_child) {
      component->child_ = Component::Create(source);
    }

    return std::move(component);
  }

  void WriteInternal(StructureSink* sink) const override {
    sink->WriteVector(vec_);
    sink->WritePoint(point_);
    sink->WriteMatrix(mat_);
    sink->WritePod<bool>(child_ != nullptr);
    if (child_) {
      child_->Write(sink);
    }
  }

  bool operator==(const Component& other) const override {
    if (!IsComponentEqualTo(other)) {
      return false;
    }
    const ComponentB& other_b = static_cast<const ComponentB&>(other);

    bool primitives_eq = (vec_ == other_b.vec_) && (point_ == other_b.point_) &&
                         (mat_ == other_b.mat_);
    if (!primitives_eq) {
      return false;
    }
    if ((child_ == nullptr) != (other_b.child_ == nullptr)) {
      return false;
    }
    if (child_ && other_b.child_) {
      if (*child_ != *other_b.child_) {
        return false;
      }
    }
    return true;
  }

  Renderable* GetRenderable() const override { return renderable_.get(); }

 private:
  std::unique_ptr<FakeRenderable> renderable_;
  Vector3i vec_;
  Point2i point_;
  Matrix4f mat_;
  std::unique_ptr<const Component> child_;

  DECLARE_SEURAT_COMPONENT(ComponentB);
};
REGISTER_SEURAT_COMPONENT(ComponentB);

// Basic label API test
TEST(ComponentTest, LabelTests) {
  ComponentB b;
  EXPECT_EQ(b.GetLabel(), std::string{});
  const std::string kNewLabel("I am b.");
  b.SetLabel(kNewLabel);
  EXPECT_EQ(b.GetLabel(), kNewLabel);
}

// Reads and writes a simple nested hierarchy of Component instances.
TEST(ComponentTest, WriteReadComponents) {
  // Create the hierarchy.
  std::unique_ptr<ComponentB> base(new ComponentB());
  {
    std::unique_ptr<ComponentA> a(new ComponentA());
    a->SetIValue(5);
    a->SetFValue(42.13f);
    std::unique_ptr<ComponentB> b(new ComponentB());
    b->SetChild(std::move(a));
    base->SetChild(std::move(b));
  }

  // Write out the Component instances.
  std::string buf;
  base::StringByteSink byte_sink(&buf);
  StructureSink sink(&byte_sink);
  base->Write(&sink);

  // Read the Component instances back out.
  size_t buffer_size = buf.size();
  base::ArrayByteSource byte_source(absl::string_view(buf.data(), buffer_size));
  StructureSource source(&byte_source);
  const std::unique_ptr<const Component> copy = Component::Create(&source);

  // Reconstruct the Component hierarchy.
  const ComponentB* const b = static_cast<const ComponentB*>(base->GetChild());
  ASSERT_TRUE(b != nullptr);
  const ComponentA* const a = static_cast<const ComponentA*>(b->GetChild());
  ASSERT_TRUE(a != nullptr);

  const ComponentB* const base_copy =
      static_cast<const ComponentB*>(copy.get());
  ASSERT_TRUE(base_copy != nullptr);
  const ComponentB* const b_copy =
      static_cast<const ComponentB*>(base_copy->GetChild());
  ASSERT_TRUE(b_copy != nullptr);
  const ComponentA* const a_copy =
      static_cast<const ComponentA*>(b_copy->GetChild());
  ASSERT_TRUE(a_copy != nullptr);

  // Validate!
  EXPECT_EQ(*base, *base_copy);
  EXPECT_EQ(*b, *b);
  EXPECT_EQ(*a, *a);
  EXPECT_EQ(base->GetVec(), base_copy->GetVec());
  EXPECT_EQ(base->GetPoint(), base_copy->GetPoint());
  EXPECT_EQ(base->GetMatrix(), base_copy->GetMatrix());
  EXPECT_EQ(b->GetVec(), b_copy->GetVec());
  EXPECT_EQ(b->GetPoint(), b_copy->GetPoint());
  EXPECT_EQ(b->GetMatrix(), b_copy->GetMatrix());
  EXPECT_EQ(a->GetIValue(), a_copy->GetIValue());
  EXPECT_EQ(a->GetFValue(), a_copy->GetFValue());
}

}  // namespace
}  // namespace component
}  // namespace seurat
