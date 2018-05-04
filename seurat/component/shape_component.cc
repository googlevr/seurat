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

#include "seurat/component/shape_component.h"

#include <limits>
#include <utility>

#include "ion/base/logging.h"
#include "ion/gfx/shaderinputregistry.h"
#include "ion/math/vector.h"
#include "seurat/component/shape_renderable.h"

namespace seurat {
namespace component {

using base::StructureSink;
using base::StructureSource;
using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Range1f;
using ion::math::Range3f;
using ion::math::Vector2f;
using ion::math::Vector2i;
using ion::math::Vector3f;

namespace {

bool ReadPrimitiveType(StructureSource* source,
                       ion::gfx::Shape::PrimitiveType* shape_type) {
  const int ion_shape_type = source->ReadPod<int>();
  if (ion_shape_type < 0 ||
      ion_shape_type > ion::gfx::Shape::PrimitiveType::kTriangleStrip) {
    return false;
  }

  *shape_type = static_cast<ion::gfx::Shape::PrimitiveType>(ion_shape_type);
  return true;
}

Range3f ExtentFromVertices(
    const std::vector<ShapeComponent::VertexAttributes>& vertices) {
  Range3f vertex_extent = std::accumulate(
      vertices.begin(), vertices.end(),
      ion::math::Range3f(vertices[0].position, vertices[0].position),
      [](const Range3f& extent,
         const ShapeComponent::VertexAttributes& vertex) {
        Range3f merged_extent = extent;
        merged_extent.ExtendByPoint(vertex.position);
        return merged_extent;
      });
  return vertex_extent;
}

}  // namespace

REGISTER_SEURAT_COMPONENT(ShapeComponent);

ShapeComponent::ShapeComponent(
    std::string label, ion::gfx::Shape::PrimitiveType ion_shape_type,
    std::vector<ShapeComponent::VertexAttributes> vertices)
    : Component(std::move(label)),
      bounding_box_(ExtentFromVertices(vertices)),
      ion_shape_type_(ion_shape_type),
      vertices_(std::move(vertices)) {}

ShapeComponent::ShapeComponent(
    std::string label, ion::gfx::Shape::PrimitiveType ion_shape_type,
    std::vector<ShapeComponent::VertexAttributes> vertices,
    std::vector<uint32> indices)
    : Component(std::move(label)),
      bounding_box_(ExtentFromVertices(vertices)),
      ion_shape_type_(ion_shape_type),
      vertices_(std::move(vertices)),
      indices_(std::move(indices)) {}

ShapeComponent::~ShapeComponent() {}

bool ShapeComponent::operator==(const Component& other) const {
  if (!Component::IsComponentEqualTo(other)) {
    return false;
  }
  const ShapeComponent& other_component =
      static_cast<const ShapeComponent&>(other);
  return (ion_shape_type_ == other_component.ion_shape_type_ &&
          vertices_ == other_component.vertices_ &&
          indices_ == other_component.indices_);
}

Renderable* ShapeComponent::GetRenderable() const { return renderable_.get(); }

Range3f ShapeComponent::GetBoundingBox() const { return bounding_box_; }

// static
std::unique_ptr<const ShapeComponent> ShapeComponent::Create(
    std::string label, ion::gfx::Shape::PrimitiveType ion_shape_type,
    std::vector<VertexAttributes> vertices) {
  std::unique_ptr<ShapeComponent> component(new ShapeComponent(
      std::move(label), ion_shape_type, std::move(vertices)));
  component->renderable_ = ShapeRenderable::Create(component.get());

  return std::move(component);
}

// static
std::unique_ptr<const ShapeComponent> ShapeComponent::Create(
    std::string label, ion::gfx::Shape::PrimitiveType ion_shape_type,
    std::vector<VertexAttributes> vertices, std::vector<uint32> indices) {
  std::unique_ptr<ShapeComponent> component(
      new ShapeComponent(std::move(label), ion_shape_type, std::move(vertices),
                         std::move(indices)));
  component->renderable_ = ShapeRenderable::Create(component.get());

  return std::move(component);
}

// static
std::unique_ptr<const Component> ShapeComponent::Create(
    std::string label, StructureSource* source) {
  // Read the Component parameters.
  ion::gfx::Shape::PrimitiveType ion_shape_type =
      ion::gfx::Shape::PrimitiveType::kTriangles;
  if (!ReadPrimitiveType(source, &ion_shape_type)) {
    return nullptr;
  }

  std::vector<VertexAttributes> vertices;
  source->ReadPodArray(&vertices);
  std::remove_const<decltype(ShapeComponent::indices_)>::type indices;
  source->ReadPodArray(&indices);

  std::unique_ptr<ShapeComponent> component(
      new ShapeComponent(std::move(label), ion_shape_type, std::move(vertices),
                         std::move(indices)));
  component->renderable_ = ShapeRenderable::Create(component.get());

  return std::move(component);
}

void ShapeComponent::WriteInternal(StructureSink* sink) const {
  sink->WritePod<int>(ion_shape_type_);
  sink->WritePodArray(vertices_);
  sink->WritePodArray(indices_);
}

}  // namespace component
}  // namespace seurat
