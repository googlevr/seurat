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

#include "seurat/mesh/mesh_renderable.h"

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
#include "seurat/component/renderable_util.h"
#include "seurat/mesh/mesh_renderable_util.h"

ION_REGISTER_ASSETS(MeshComponentAssets);

namespace seurat {
namespace mesh {

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
    const MeshComponent& mesh_component,
    const ion::gfx::ShaderInputRegistryPtr& registry) {
  // Create and fill attribute buffer and index buffer.
  ion::gfx::BufferObjectPtr attribute_buffer =
      component::CreateAttributeBufferFromVector(
          mesh_component.GetAttributeBuffer());
  ion::gfx::IndexBufferPtr index_buffer =
      component::CreateIndexBufferFromVector(
          mesh_component.GetAttributeBuffer().size(),
          mesh_component.GetIndexBuffer());

  // Bind the attribute buffer to an attribute array.
  ion::gfx::AttributeArrayPtr attribute_array(new ion::gfx::AttributeArray);
  MeshComponent::VertexAttributes sa;
  ion::gfxutils::BufferToAttributeBinder<MeshComponent::VertexAttributes>(sa)
      .Bind(sa.position, "aPosition")
      .Bind(sa.tex_coord, "aTexCoord")
      .Apply(registry, attribute_array, attribute_buffer);

  // Create the shape and add it to the node.
  ion::gfx::ShapePtr shape(new ion::gfx::Shape);
  shape->SetPrimitiveType(ion::gfx::Shape::kTriangles);
  shape->SetAttributeArray(attribute_array);
  shape->SetIndexBuffer(index_buffer);

  return shape;
}

ion::gfx::NodePtr CreateMeshNode(const Renderable::RenderingContext& context,
                                 const MeshComponent& component,
                                 size_t* u_clip_from_mesh_matrix_index) {
  struct MeshShaderSpec {
    std::string name;
    std::string vertex_shader_path;
    std::string fragment_shader_path;
  };

  const MeshShaderSpec kTexturedMeshShader{"TexturedMesh", "shaders/mesh.vert",
                                           "shaders/mesh.frag"};
  const MeshShaderSpec kSolidDepthMeshShader{"SolidDepthMesh",
                                             "shaders/solid_depth_mesh.vert",
                                             "shaders/solid_depth_mesh.frag"};

  MeshShaderSpec shader_spec;
  // Switch shader depending on whether there is a texture.
  if (component.HasTextureAtlas()) {
    shader_spec = kTexturedMeshShader;
  } else {
    shader_spec = kSolidDepthMeshShader;
  }

  // Either reuse the shader registry and program of an existing opaque shader
  // program (created by another mesh renderable), or create a new one.
  ion::gfx::ShaderInputRegistryPtr shader_registry;
  ion::gfx::ShaderProgramPtr shader_program;
  component::LazyConstructShaderRegistryAndProgram(
      context.shader_manager, shader_spec.name, shader_spec.vertex_shader_path,
      shader_spec.fragment_shader_path, &shader_registry, &shader_program);

  ion::gfx::ShapePtr shape = CreateShape(component, shader_registry);

  ion::gfx::NodePtr node =
      BuildRenderNode(context, shader_registry, shader_program, shape,
                      u_clip_from_mesh_matrix_index);

  // If a texture exists, bind it to uTexture.
  if (component.HasTextureAtlas()) {
    ion::gfx::TexturePtr color_texture =
        base::CreateTexture(component.GetTextureAtlas());
    // Use bilinear texture filtering.
    auto color_texture_sampler = color_texture->GetSampler();
    color_texture_sampler->SetMinFilter(ion::gfx::Sampler::FilterMode::kLinear);
    color_texture_sampler->SetMagFilter(ion::gfx::Sampler::FilterMode::kLinear);
    node->AddUniform(
        shader_registry->Create<ion::gfx::Uniform>("uTexture", color_texture));
  }
  return node;
}

}  // namespace

MeshRenderable::MeshRenderable(const MeshComponent* mesh_component)
    : mesh_component_(mesh_component) {
  ION_STATIC_ONCE_CHECKED(MeshComponentAssets::RegisterAssets);
}

Renderable::Properties MeshRenderable::GetProperties() const {
  Properties properties;
  properties.bounding_box = mesh_component_->GetBoundingBox();
  return properties;
}

Renderable::NodeSet MeshRenderable::GetRenderNodes(
    const RenderingContext& context) {
  // Return an empty node set if the component is empty.
  if (mesh_component_->IsEmpty()) {
    NodeSet node_set;
    node_set.opaque.Reset(new ion::gfx::Node);
    node_set.transparent.Reset(new ion::gfx::Node);
    return node_set;
  }

  // Initialize node_.
  if (node_.Get() == nullptr) {
    node_ = CreateMeshNode(context, *mesh_component_,
                           &u_clip_from_mesh_matrix_index_);
  }

  NodeSet node_set;

  // Opaque
  {
    node_set.opaque.Reset(new ion::gfx::Node);
    // If there is no texture, render to the depth-buffer only.
    if (!mesh_component_->HasTextureAtlas()) {
      node_set.opaque->AddChild(node_);

      ion::gfx::StateTablePtr state(new ion::gfx::StateTable);
      state->Enable(ion::gfx::StateTable::Capability::kStencilTest, false);
      state->Enable(ion::gfx::StateTable::Capability::kCullFace, false);
      state->Enable(ion::gfx::StateTable::Capability::kBlend, false);
      state->SetDepthWriteMask(true);
      state->Enable(ion::gfx::StateTable::kDepthTest, true);
      state->Enable(ion::gfx::StateTable::kSampleAlphaToCoverage, false);
      state->SetColorWriteMasks(false, false, false, false);
      node_set.opaque->SetStateTable(state);
    }
  }

  // Transparent
  {
    node_set.transparent.Reset(new ion::gfx::Node);
    // If there is an alpha texture, render to color-only with
    // premultiplied-alpha blending.
    if (mesh_component_->HasTextureAtlas()) {
      const bool kIsAlphaPremultiplied = true;
      node_set.transparent =
          BuildTransparentPassNode(node_, kIsAlphaPremultiplied);
    }
  }

  return node_set;
}

void MeshRenderable::Update(const ViewState& view_state) {
  // Early out if the component is empty.
  if (mesh_component_->IsEmpty()) return;

  CHECK_NOTNULL(node_.Get());
  node_->SetUniformValueAt(u_clip_from_mesh_matrix_index_, 0,
                           view_state.left_clip_from_component_matrix);
  node_->SetUniformValueAt(u_clip_from_mesh_matrix_index_, 1,
                           view_state.right_clip_from_component_matrix);
}

}  // namespace mesh
}  // namespace seurat
