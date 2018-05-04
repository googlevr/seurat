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

#include "seurat/component/shape_component_utils.h"

namespace seurat {
namespace component {

using base::Color4f;
using ion::math::Point3f;

std::vector<ShapeComponent::VertexAttributes> MakeWireframeBoxVertices(
    const Point3f& minimal, const Point3f& maximal, const Color4f& color) {
  std::vector<component::ShapeComponent::VertexAttributes> box_edges;

  for (int axis0 = 0; axis0 < 3; ++axis0) {
    // Run this algorithm over each axis:
    // Iterate the four corners of the face, and construct the edge at each
    // corner, spanning minimal to maximal range of the box on that axis.
    int axis1 = (axis0 + 1) % 3;
    int axis2 = (axis1 + 1) % 3;
    // Loop the corners of the face on the axis1-axis2 face.
    for (int axis1_end = 0; axis1_end < 2; ++axis1_end) {
      for (int axis2_end = 0; axis2_end < 2; ++axis2_end) {
        Point3f endpoint;
        // Start at the minimal endpoint on the outer axis.
        endpoint[axis0] = minimal[axis0];
        // Fill in the corner values.
        endpoint[axis1] = (axis1_end ? maximal : minimal)[axis1];
        endpoint[axis2] = (axis2_end ? maximal : minimal)[axis2];
        box_edges.push_back({endpoint, color});

        // Span to the maximal endpoint on the outer axis.
        endpoint[axis0] = maximal[axis0];
        box_edges.push_back({endpoint, color});
      }
    }
  }

  return box_edges;
}

std::unique_ptr<const ShapeComponent> MakeWireframeBoxShape(
    std::string label, const ion::math::Point3f& minimal,
    const ion::math::Point3f& maximal, const base::Color4f& color) {
  return ShapeComponent::Create(
      std::move(label), ion::gfx::Shape::PrimitiveType::kLines,
      MakeWireframeBoxVertices(minimal, maximal, color));
}

}  // namespace component
}  // namespace seurat
