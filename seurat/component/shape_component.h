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

#ifndef VR_SEURAT_COMPONENT_SHAPE_COMPONENT_H_
#define VR_SEURAT_COMPONENT_SHAPE_COMPONENT_H_

#include <memory>
#include <vector>

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "seurat/base/structured_io.h"
#include "seurat/base/util.h"
#include "seurat/component/component.h"

namespace seurat {
namespace component {

// ShapeComponent supports serializing and rendering a general Ion Shape,
// intended for visualizing data.
class ShapeComponent : public Component {
 public:
  // Defines a single vertex of the shape.
  struct VertexAttributes {
    ion::math::Point3f position;
    base::Color4f color;

    bool operator==(const VertexAttributes& rhs) const {
      return rhs.position == position && rhs.color == color;
    }

    bool operator!=(const VertexAttributes& rhs) const {
      return !operator==(rhs);
    }
  };

  // Constructs a shape without using vertex indices.
  ShapeComponent(std::string label,
                 ion::gfx::Shape::PrimitiveType ion_shape_type,
                 std::vector<VertexAttributes> vertices);

  // Constructs a shape from |vertices| indexed by |indices|. Labels the Ion
  // node with human-readable description |label|.
  ShapeComponent(std::string label,
                 ion::gfx::Shape::PrimitiveType ion_shape_type,
                 std::vector<VertexAttributes> vertices,
                 std::vector<uint32> indices);

  ~ShapeComponent() override;

  // Creates a shape component for rendering, of topology |ion_shape_type| and
  // vertex data |vertices|, but no indices. Labels the Ion node with optional
  // human-readable description |label|.
  static std::unique_ptr<const ShapeComponent> Create(
      std::string label, ion::gfx::Shape::PrimitiveType ion_shape_type,
      std::vector<VertexAttributes> vertices);

  // Creates a shape component for rendering, of topology |ion_shape_type| and
  // vertex data |vertices|, with vertices indexed by |indices|. Labels the Ion
  // node with optional human-readable description |label|.
  static std::unique_ptr<const ShapeComponent> Create(
      std::string label, ion::gfx::Shape::PrimitiveType ion_shape_type,
      std::vector<VertexAttributes> vertices, std::vector<uint32> indices);

  static std::unique_ptr<const Component> Create(std::string label,
                                                 base::StructureSource* source);

  // Component implementation.
  bool operator==(const Component& other) const override;
  Renderable* GetRenderable() const override;

  // Returns true if this component has no vertices.
  bool IsEmpty() const { return vertices_.empty(); }

  // Calculates the extent of the entire component.
  ion::math::Range3f GetBoundingBox() const;

  // Returns the shape's vertices.
  const std::vector<VertexAttributes>& GetVertices() const { return vertices_; }

  // Returns the index buffer, which may be empty.
  const std::vector<uint32>& GetIndices() const { return indices_; }

  // Returns the primitive topology.
  ion::gfx::Shape::PrimitiveType GetPrimitiveType() const {
    return ion_shape_type_;
  }

 private:
  // Component implementation.
  void WriteInternal(base::StructureSink* sink) const override;

  // The 3D extent of the entire component.
  const ion::math::Range3f bounding_box_;

  // The shape's primitive topology, i.e. the interpretation of vertex
  // connectivity.
  ion::gfx::Shape::PrimitiveType ion_shape_type_;

  // Stores the vertices of the shape.
  const std::vector<VertexAttributes> vertices_;

  // Stores optional indices describing connectivity of the shape; otherwise
  // vertex order defines connectivity. See ion::gfx::Shape::SetIndexBuffer.
  const std::vector<uint32> indices_;

  // The Renderable which renders the contents of this Component.
  std::unique_ptr<Renderable> renderable_;

  DECLARE_SEURAT_COMPONENT(ShapeComponent);
  DISALLOW_COPY_AND_ASSIGN(ShapeComponent);
};

// Enables VertexAttributes bitwise serialization.
inline void WriteObjectsToSink(const ShapeComponent::VertexAttributes* objects,
                               size_t count, base::StructureSink* sink) {
  sink->WriteObjects(objects, count);
}

// Enables VertexAttributes bitwise deserialization.
inline void ReadObjectsFromSource(base::StructureSource* source,
                                  ShapeComponent::VertexAttributes* objects,
                                  size_t count) {
  source->ReadObjects(objects, count);
}

}  // namespace component
}  // namespace seurat

#endif  // VR_SEURAT_COMPONENT_SHAPE_COMPONENT_H_
