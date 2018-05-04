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

#include "seurat/component/shape_renderable.h"

#include <memory>
#include <string>

#include "ion/base/allocator.h"
#include "ion/base/datacontainer.h"
#include "ion/base/once.h"
#include "ion/base/zipassetmanagermacros.h"
#include "ion/gfx/attributearray.h"
#include "ion/gfx/bufferobject.h"
#include "ion/gfx/graphicsmanager.h"
#include "ion/gfx/image.h"
#include "ion/gfx/node.h"
#include "ion/gfx/shaderinputregistry.h"
#include "ion/gfx/shape.h"
#include "ion/gfx/texture.h"
#include "ion/gfx/uniform.h"
#include "ion/gfxutils/buffertoattributebinder.h"
#include "ion/gfxutils/shadermanager.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/component/renderable_util.h"
#include "seurat/component/shape_component.h"

ION_REGISTER_ASSETS(ShapeAssets);

namespace seurat {
namespace component {

using ion::math::Matrix4f;

namespace {

const char kShapeShaderName[] = "Shape";
const char kShapeVertexShaderName[] = "shaders/shape.vert";
const char kShapeFragmentShaderName[] = "shaders/shape.frag";

}  // namespace

ShapeRenderable::ShapeRenderable(const ShapeComponent* shape_component)
    : shape_component_(shape_component) {
  ShapeAssets::RegisterAssetsOnce();
}

ShapeRenderable::~ShapeRenderable() {}

// static
std::unique_ptr<ShapeRenderable> ShapeRenderable::Create(
    const ShapeComponent* shape_component) {
  return std::unique_ptr<ShapeRenderable>(new ShapeRenderable(shape_component));
}

Renderable::Properties ShapeRenderable::GetProperties() const {
  Properties properties;
  properties.bounding_box = shape_component_->GetBoundingBox();
  return properties;
}

Renderable::NodeSet ShapeRenderable::GetRenderNodes(
    const RenderingContext& context) {
  // Return an empty node set if the component is empty.
  if (shape_component_->IsEmpty()) {
    node_set_.transparent.Reset(new ion::gfx::Node);
    return node_set_;
  }

  // The shader registry is shared between the soft-z and color accumulation
  // shader programs, because they draw the same shape. We either reuse the
  // registry of an existing soft-z shader program (created by another shape
  // renderable), or we create a new one.
  ion::gfx::ShaderInputRegistryPtr shader_registry;
  ion::gfxutils::ShaderManagerPtr shader_manager = context.shader_manager;
  ion::gfx::ShaderProgramPtr utility_shader_program;
  LazyConstructShaderRegistryAndProgram(
      shader_manager, kShapeShaderName, kShapeVertexShaderName,
      kShapeFragmentShaderName, &shader_registry, &utility_shader_program);

  if (shape_node_.Get() == nullptr) {
    shape_node_.Reset(new ion::gfx::Node);
    shape_node_->SetLabel(shape_component_->GetLabel());
  }
  if (shape_.Get() == nullptr) {
    shape_.Reset(new ion::gfx::Shape);
    shape_->SetPrimitiveType(shape_component_->GetPrimitiveType());

    shape_node_->SetShaderProgram(utility_shader_program);
    if (node_set_.transparent.Get() == nullptr) {
      node_set_.transparent.Reset(new ion::gfx::Node);
    }
    node_set_.transparent->AddChild(shape_node_);

    // Create and fill attribute buffer and index buffer.
    ion::gfx::BufferObjectPtr attribute_buffer =
        CreateAttributeBufferFromVector(shape_component_->GetVertices());

    ion::gfx::IndexBufferPtr index_buffer;
    if (!shape_component_->GetIndices().empty()) {
      index_buffer =
          CreateIndexBufferFromVector(shape_component_->GetVertices().size(),
                                      shape_component_->GetIndices());
    }

    // Bind the attribute buffer to an attribute array.
    ion::gfx::AttributeArrayPtr attribute_array(new ion::gfx::AttributeArray);
    ShapeComponent::VertexAttributes sa;
    ion::gfxutils::BufferToAttributeBinder<ShapeComponent::VertexAttributes>(sa)
        .Bind(sa.position, "aPosition")
        .Bind(sa.color, "aColor")
        .Apply(shader_registry, attribute_array, attribute_buffer);

    // Create the shape and add it to the node.
    shape_->SetAttributeArray(attribute_array);
    shape_->SetIndexBuffer(index_buffer);

    shape_node_->AddShape(shape_);
  }

  CreateUniforms(shader_registry);

  return node_set_;
}

void ShapeRenderable::CreateUniforms(
    const ion::gfx::ShaderInputRegistryPtr& shader_registry) {
  // Create uniforms for the clip-from-shape space matrices.  The matrices are
  // uninitialized, because we don't know the clip-from-component matrices yet.
  // Correct values will be set by Update() for every frame.

  u_leftright_clip_from_shape_matrix_index_ =
      shape_node_->AddUniform(shader_registry->CreateArrayUniform<Matrix4f>(
          "uClipFromShapeMatrix", nullptr, 2, ion::base::AllocatorPtr()));
  CHECK_NE(u_leftright_clip_from_shape_matrix_index_, ion::base::kInvalidIndex);
}

void ShapeRenderable::Update(const ViewState& view_state) {
  // Early out if the component is empty.
  if (shape_.Get() == nullptr) return;

  CHECK_NOTNULL(shape_node_.Get());
  shape_node_->SetUniformValueAt(u_leftright_clip_from_shape_matrix_index_, 0,
                                 view_state.left_clip_from_component_matrix);
  shape_node_->SetUniformValueAt(u_leftright_clip_from_shape_matrix_index_, 1,
                                 view_state.right_clip_from_component_matrix);
}

}  // namespace component
}  // namespace seurat
