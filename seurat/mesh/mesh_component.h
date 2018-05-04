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

#ifndef VR_SEURAT_MESH_MESH_COMPONENT_H_
#define VR_SEURAT_MESH_MESH_COMPONENT_H_

#include <array>

#include "ion/gfx/image.h"
#include "ion/math/vector.h"
#include "seurat/component/component.h"

namespace seurat {
namespace mesh {

// A component for rendering a textured triangle mesh.
class MeshComponent : public component::Component {
 public:
  struct VertexAttributes {
    VertexAttributes() = default;
    VertexAttributes(const ion::math::Point3f& position,
                     const ion::math::Point3f& tex_coord)
        : position(position), tex_coord(tex_coord) {}

    bool operator==(const VertexAttributes& rhs) const {
      return position == rhs.position && tex_coord == rhs.tex_coord;
    }

    ion::math::Point3f position;
    ion::math::Point3f tex_coord;
  };

  ~MeshComponent() override = default;

  // Creates a new MeshComponent with the given texture atlas, attribute buffer,
  // and index buffer.
  static std::unique_ptr<MeshComponent> Create(
      std::string label, ion::gfx::ImagePtr texture_atlas,
      std::vector<VertexAttributes> attribute_buf,
      std::vector<uint32> index_buf);

  // Component implementation.
  static std::unique_ptr<const component::Component> Create(
      std::string label, base::StructureSource* source);
  component::Renderable* GetRenderable() const override;
  bool operator==(const component::Component& other) const override;

  // Returns true if this mesh component contains no geometry.
  bool IsEmpty() const { return attribute_buf_.empty() || index_buf_.empty(); }

  // Returns true if this mesh component has a texture atlas.
  bool HasTextureAtlas() const { return GetTextureAtlas().Get() != nullptr; }

  // Returns the texture for this mesh.
  const ion::gfx::ImagePtr& GetTextureAtlas() const { return texture_atlas_; }

  // Returns the attribute buffer for this mesh, suitable for sending directly
  // to OpenGL.
  const std::vector<VertexAttributes>& GetAttributeBuffer() const {
    return attribute_buf_;
  }

  // Returns the index buffer, suitable for sending directly to OpenGL.
  const std::vector<uint32>& GetIndexBuffer() const { return index_buf_; }

  // Returns an axis-aligned bounding box for the triangles in the mesh.
  const ion::math::Range3f& GetBoundingBox() const { return aabb_; }

 private:
  MeshComponent(std::string label, ion::gfx::ImagePtr texture_atlas,
                std::vector<VertexAttributes> attribute_buf,
                std::vector<uint32> index_buf, const ion::math::Range3f& aabb);

  // Component implementation.
  void WriteInternal(base::StructureSink* sink) const override;

  // The texture atlas.
  const ion::gfx::ImagePtr texture_atlas_;

  // Interleaved attribute-buffer which can be sent directly to OpenGL for
  // rendering.
  const std::vector<VertexAttributes> attribute_buf_;

  // The index buffer associated with attribute_buf_ for use with GL_TRIANGLES.
  //
  // Unsigned integers are used for direct compatibility with OpenGL.
  const std::vector<uint32> index_buf_;

  // Axis-aligned bounding box for the mesh.
  const ion::math::Range3f aabb_;

  // The renderable used for rendering this component.
  const std::unique_ptr<component::Renderable> renderable_;

  DECLARE_SEURAT_COMPONENT(MeshComponent);
};

}  // namespace mesh
}  // namespace seurat

#endif  // VR_SEURAT_MESH_MESH_COMPONENT_H_
