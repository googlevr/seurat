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

#ifndef VR_SEURAT_COMPONENT_RENDERABLE_UTIL_H_
#define VR_SEURAT_COMPONENT_RENDERABLE_UTIL_H_

#include <limits>
#include <string>
#include <vector>

#include "ion/base/datacontainer.h"
#include "ion/gfx/bufferobject.h"
#include "ion/gfx/indexbuffer.h"
#include "ion/gfx/shaderinputregistry.h"
#include "ion/gfx/shaderprogram.h"
#include "ion/gfxutils/buffertoattributebinder.h"
#include "ion/gfxutils/shadermanager.h"
#include "absl/types/span.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/component/renderable.h"

namespace seurat {
namespace component {

// Creates an Ion DataContainer from a std::vector compatible-container of
// arbitrary value type.
template <typename Container>
ion::base::DataContainerPtr DataContainerFromVector(Container data_vector) {
  return ion::base::DataContainer::CreateAndCopy<
      typename Container::value_type>(data_vector.data(), data_vector.size(),
                                      false, ion::base::AllocatorPtr());
}

// Creates and returns an Ion BufferObject from the given |attribute_vector|.
template <typename Container>
ion::gfx::BufferObjectPtr CreateAttributeBufferFromVector(
    Container attribute_vector) {
  ion::gfx::BufferObjectPtr attribute_buffer(new ion::gfx::BufferObject);
  // Create and populate the attribute buffer.
  ion::base::DataContainerPtr attribute_container =
      DataContainerFromVector(attribute_vector);
  attribute_buffer->SetData(
      attribute_container, sizeof(typename Container::value_type),
      attribute_vector.size(), ion::gfx::BufferObject::kStaticDraw);
  return attribute_buffer;
}

// Creates and returns an Ion IndexBuffer from the given |index_vector|. The
// number of vertices must be specified, because the implementation uses 16-bit
// or 32-bit indices, depending on the number of vertices.
ion::gfx::IndexBufferPtr CreateIndexBufferFromVector(
    int num_vertices, absl::Span<const uint32> index_vector);

// Lazily constructs the |shader_registry| and creates the |shader_program| if
// it doesn't already exist.
void LazyConstructShaderRegistryAndProgram(
    const ion::gfxutils::ShaderManagerPtr& shader_manager,
    const std::string& shader_name, const std::string& vertex_shader_path,
    const std::string& fragment_shader_path,
    ion::gfx::ShaderInputRegistryPtr* shader_registry,
    ion::gfx::ShaderProgramPtr* shader_program);

// Helps construct Ion Shapes.
//
// Only works with GL_TRIANGLES & non-interleaved attributes.
//
// Usage:
//   ShapeBuilder builder(registry);
//   builder.BindAttributeFromVector("uPosition", positions);
//   builder.BindAttributeFromVector("uTexCoord", tex_coords);
//   shape = builder.Build(positions.size(), indices);
class ShapeBuilder {
 public:
  explicit ShapeBuilder(const ion::gfx::ShaderInputRegistryPtr& registry)
      : registry_(registry), attribute_array_(new ion::gfx::AttributeArray) {}

  // Binds a vector<AttributeT> to the attribute with the specified |name|,
  // copying the data into an ion buffer.
  template <typename AttributeT>
  void BindAttributeFromVector(
      const std::string& name,
      const std::vector<AttributeT>& attribute_vector) {
    AttributeT a;
    ion::gfxutils::BufferToAttributeBinder<AttributeT>(a)  //
        .Bind(a, name)                                     //
        .Apply(registry_, attribute_array_,
               CreateAttributeBufferFromVector(attribute_vector));
  }

  // Builds a Shape including all previously-bound attributes.
  ion::gfx::ShapePtr Build(int vertex_count,
                           const std::vector<uint32>& index_buffer);

 private:
  const ion::gfx::ShaderInputRegistryPtr registry_;
  const ion::gfx::AttributeArrayPtr attribute_array_;
};

}  // namespace component
}  // namespace seurat

#endif  // VR_SEURAT_COMPONENT_RENDERABLE_UTIL_H_
