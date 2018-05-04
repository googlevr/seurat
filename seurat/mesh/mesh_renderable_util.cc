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

#include "seurat/mesh/mesh_renderable_util.h"

#include <array>

#include "ion/math/matrix.h"

namespace seurat {
namespace mesh {

using component::Renderable;
using ion::math::Matrix4f;

ion::gfx::NodePtr BuildRenderNode(
    const Renderable::RenderingContext &context,
    const ion::gfx::ShaderInputRegistryPtr &shader_registry,
    const ion::gfx::ShaderProgramPtr &shader_program,
    const ion::gfx::ShapePtr &shape, size_t *u_clip_from_mesh_matrix_index) {
  // A node which will render our mesh.
  ion::gfx::NodePtr node(new ion::gfx::Node);
  node->SetShaderProgram(shader_program);
  node->AddShape(shape);

  std::array<Matrix4f, 2> initial_values{
      {Matrix4f::Identity(), Matrix4f::Identity()}};
  *u_clip_from_mesh_matrix_index =
      node->AddUniform(shader_registry->CreateArrayUniform<Matrix4f>(
          "uClipFromMeshMatrix", initial_values.data(), 2,
          ion::base::AllocatorPtr()));
  CHECK_NE(*u_clip_from_mesh_matrix_index, ion::base::kInvalidIndex);

  return node;
}

ion::gfx::NodePtr BuildTransparentPassNode(const ion::gfx::NodePtr &node,
                                           bool is_alpha_premultiplied) {
  ion::gfx::NodePtr transparent(new ion::gfx::Node);
  transparent->AddChild(node);
  ion::gfx::StateTablePtr state(new ion::gfx::StateTable);
  state->Enable(ion::gfx::StateTable::Capability::kCullFace, false);
  state->Enable(ion::gfx::StateTable::Capability::kBlend, true);
  state->SetDepthWriteMask(false);
  state->Enable(ion::gfx::StateTable::kDepthTest, true);
  state->SetBlendEquations(ion::gfx::StateTable::BlendEquation::kAdd,
                           ion::gfx::StateTable::BlendEquation::kAdd);
  ion::gfx::StateTable::BlendFunctionFactor rgb_src_factor =
      is_alpha_premultiplied
          ? ion::gfx::StateTable::BlendFunctionFactor::kOne
          : ion::gfx::StateTable::BlendFunctionFactor::kSrcAlpha;
  state->SetBlendFunctions(
      rgb_src_factor,
      ion::gfx::StateTable::BlendFunctionFactor::kOneMinusSrcAlpha,
      ion::gfx::StateTable::BlendFunctionFactor::kOne,
      ion::gfx::StateTable::BlendFunctionFactor::kZero);
  transparent->SetStateTable(state);
  return transparent;
}

}  // namespace mesh
}  // namespace seurat
