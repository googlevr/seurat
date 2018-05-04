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

#include "seurat/component/renderable_util.h"

#include <array>

#include "ion/math/matrix.h"
#include "absl/types/span.h"

namespace seurat {

namespace component {

ion::gfx::IndexBufferPtr CreateIndexBufferFromVector(
    int num_vertices, absl::Span<const uint32> index_vector) {
  ion::gfx::IndexBufferPtr index_buffer(new ion::gfx::IndexBuffer);
  // Create and populate the index-buffer.
  // Use 16-bit indices if possible.
  if (num_vertices < std::numeric_limits<uint16>::max()) {
    // Convert from uint32 to uint16.
    std::vector<uint16> index_vector_uint16(index_vector.begin(),
                                            index_vector.end());

    ion::base::DataContainerPtr index_container =
        DataContainerFromVector(index_vector_uint16);
    index_buffer->AddSpec(ion::gfx::BufferObject::kUnsignedShort, 1, 0);
    index_buffer->SetData(index_container, sizeof(uint16), index_vector.size(),
                          ion::gfx::BufferObject::kStaticDraw);
  } else {
    ion::base::DataContainerPtr index_container =
        DataContainerFromVector(index_vector);
    index_buffer->AddSpec(ion::gfx::BufferObject::kUnsignedInt, 1, 0);
    index_buffer->SetData(index_container, sizeof(uint32), index_vector.size(),
                          ion::gfx::BufferObject::kStaticDraw);
  }
  return index_buffer;
}

void LazyConstructShaderRegistryAndProgram(
    const ion::gfxutils::ShaderManagerPtr& shader_manager,
    const std::string& shader_name, const std::string& vertex_shader_path,
    const std::string& fragment_shader_path,
    ion::gfx::ShaderInputRegistryPtr* shader_registry,
    ion::gfx::ShaderProgramPtr* shader_program) {
  *shader_program = shader_manager->GetShaderProgram(shader_name);
  if ((*shader_program).Get() == nullptr) {
    (*shader_registry).Reset(new ion::gfx::ShaderInputRegistry);
    (*shader_registry)->IncludeGlobalRegistry();
    *shader_program =
        base::CreateShaderProgram(shader_manager, *shader_registry, shader_name,
                                  vertex_shader_path, fragment_shader_path);
  } else {
    *shader_registry =
        shader_manager->GetShaderProgram(shader_name)->GetRegistry();
  }
  CHECK_NOTNULL((*shader_registry).Get());
  CHECK_NOTNULL((*shader_program).Get());
}

ion::gfx::ShapePtr ShapeBuilder::Build(
    int vertex_count, const std::vector<uint32>& index_buffer) {
  ion::gfx::ShapePtr shape(new ion::gfx::Shape);
  shape->SetPrimitiveType(ion::gfx::Shape::kTriangles);
  shape->SetAttributeArray(attribute_array_);
  shape->SetIndexBuffer(
      CreateIndexBufferFromVector(vertex_count, index_buffer));
  return shape;
}

}  // namespace component
}  // namespace seurat
