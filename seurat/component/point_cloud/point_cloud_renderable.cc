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

#include "seurat/component/point_cloud/point_cloud_renderable.h"

#include <array>
#include <memory>
#include <string>

#include "ion/base/allocator.h"
#include "ion/base/once.h"
#include "ion/base/zipassetmanagermacros.h"
#include "ion/gfx/attributearray.h"
#include "ion/gfx/bufferobject.h"
#include "ion/gfx/image.h"
#include "ion/gfx/indexbuffer.h"
#include "ion/gfx/node.h"
#include "ion/gfx/shaderinputregistry.h"
#include "ion/gfx/shape.h"
#include "ion/gfx/statetable.h"
#include "ion/gfx/texture.h"
#include "ion/gfx/uniform.h"
#include "ion/gfxutils/buffertoattributebinder.h"
#include "ion/math/vector.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/component/point_cloud/point_cloud_component.h"
#include "seurat/component/renderable_util.h"

ION_REGISTER_ASSETS(PointCloudAssets);

namespace seurat {
namespace component {

using component::Renderable;
using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Point2i;
using ion::math::Point2ui16;
using ion::math::Point3f;
using ion::math::Point4ui16;
using ion::math::Vector2i;
using ion::math::Vector2ui16;
using ion::math::Vector3f;
using ion::math::Vector4ui16;

namespace {

ion::gfx::ShapePtr CreateShape(
    const PointCloudComponent& point_cloud_component,
    const ion::gfx::ShaderInputRegistryPtr& registry) {
  // Create and fill attribute buffer and index buffer.
  ion::gfx::BufferObjectPtr attribute_buffer =
      component::CreateAttributeBufferFromVector(
          point_cloud_component.GetAttributeBuffer());

  // Bind the attribute buffer to an attribute array.
  ion::gfx::AttributeArrayPtr attribute_array(new ion::gfx::AttributeArray);
  PointCloudComponent::VertexAttributes sa;
  ion::gfxutils::BufferToAttributeBinder<PointCloudComponent::VertexAttributes>(
      sa)
      .Bind(sa.position, "aPosition")
      .Apply(registry, attribute_array, attribute_buffer);

  // Create the shape and add it to the node.
  ion::gfx::ShapePtr shape(new ion::gfx::Shape);
  shape->SetPrimitiveType(ion::gfx::Shape::kPoints);
  shape->SetAttributeArray(attribute_array);

  return shape;
}

ion::gfx::NodePtr CreatePointCloudNode(
    const Renderable::RenderingContext& context,
    const PointCloudComponent& component,
    size_t* u_clip_from_object_matrix_index) {
  std::string vertex_shader_path = "point_cloud.vert";
  std::string fragment_shader_path = "point_cloud.frag";

  // Either reuse the shader registry and program of an existing opaque shader
  // program (created by another point cloud renderable), or create a new one.
  ion::gfx::ShaderInputRegistryPtr shader_registry;
  ion::gfx::ShaderProgramPtr shader_program;
  component::LazyConstructShaderRegistryAndProgram(
      context.shader_manager, "PointCloud", vertex_shader_path,
      fragment_shader_path, &shader_registry, &shader_program);

  ion::gfx::ShapePtr shape = CreateShape(component, shader_registry);

  // A node which will render our mesh.
  ion::gfx::NodePtr node(new ion::gfx::Node);
  node->SetShaderProgram(shader_program);
  node->AddShape(shape);

  std::array<Matrix4f, 2> initial_values{
      {Matrix4f::Identity(), Matrix4f::Identity()}};
  *u_clip_from_object_matrix_index =
      node->AddUniform(shader_registry->CreateArrayUniform<Matrix4f>(
          "uClipFromObjectMatrix", initial_values.data(), 2,
          ion::base::AllocatorPtr()));
  CHECK_NE(*u_clip_from_object_matrix_index, ion::base::kInvalidIndex);

  return node;
}

}  // namespace

PointCloudRenderable::PointCloudRenderable(
    const PointCloudComponent* point_cloud_component)
    : point_cloud_component_(point_cloud_component) {
  ION_STATIC_ONCE_CHECKED(PointCloudAssets::RegisterAssets);
}

Renderable::Properties PointCloudRenderable::GetProperties() const {
  Properties properties;
  properties.bounding_box = point_cloud_component_->GetBoundingBox();
  return properties;
}

Renderable::NodeSet PointCloudRenderable::GetRenderNodes(
    const RenderingContext& context) {
  NodeSet node_set;
  node_set.opaque.Reset(new ion::gfx::Node);

  // Return an empty node set if the component is empty.
  if (point_cloud_component_->IsEmpty()) {
    return node_set;
  }

  // Initialize node_.
  if (node_.Get() == nullptr) {
    node_ = CreatePointCloudNode(context, *point_cloud_component_,
                                 &u_clip_from_object_matrix_index_);
  }

  // Transparent
  {
    node_set.transparent.Reset(new ion::gfx::Node);
    node_set.transparent->AddChild(node_);
    ion::gfx::StateTablePtr state(new ion::gfx::StateTable);
    state->Enable(ion::gfx::StateTable::Capability::kStencilTest, false);
    state->Enable(ion::gfx::StateTable::Capability::kCullFace, false);
    state->Enable(ion::gfx::StateTable::Capability::kBlend, true);
    state->SetDepthWriteMask(false);
    state->SetBlendEquations(ion::gfx::StateTable::BlendEquation::kAdd,
                             ion::gfx::StateTable::BlendEquation::kAdd);
    state->SetBlendFunctions(
        ion::gfx::StateTable::BlendFunctionFactor::kOne,
        ion::gfx::StateTable::BlendFunctionFactor::kOneMinusSrcAlpha,
        ion::gfx::StateTable::BlendFunctionFactor::kOne,
        ion::gfx::StateTable::BlendFunctionFactor::kZero);

    node_set.transparent->SetStateTable(state);
  }

  return node_set;
}

void PointCloudRenderable::Update(const ViewState& view_state) {
  // Early out if the component is empty.
  if (point_cloud_component_->IsEmpty()) return;

  CHECK_NOTNULL(node_.Get());
  node_->SetUniformValueAt(u_clip_from_object_matrix_index_, 0,
                           view_state.left_clip_from_component_matrix);
  node_->SetUniformValueAt(u_clip_from_object_matrix_index_, 1,
                           view_state.right_clip_from_component_matrix);
}

}  // namespace component
}  // namespace seurat
