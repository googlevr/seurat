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

#include <memory>
#include <vector>

#include "gtest/gtest.h"
#include "seurat/base/color.h"
#include "seurat/component/component_tests_utils.h"
#include "seurat/component/shape_component.h"

namespace seurat {
namespace component {
namespace {

using base::Color4f;
using ion::math::Point3f;
using ion::math::Range3f;

TEST(ShapeComponentUtilsTest, BoxProperties) {
  const Point3f kMin(-1.0f, 2.0f, 0.0f);
  const Point3f kMax(-0.5f, 4.0f, 3.1f);
  const Color4f kYellow(1.0f, 1.0f, 0.0f, 1.0f);
  std::vector<ShapeComponent::VertexAttributes> box_edges =
      MakeWireframeBoxVertices(kMin, kMax, kYellow);

  // Check we got 2 vertices/edge * 12 edges = 24 vertices.
  EXPECT_EQ(box_edges.size(), 24);
  Range3f box = std::accumulate(
      box_edges.begin(), box_edges.end(), Range3f(),
      [](const Range3f& lhs, const ShapeComponent::VertexAttributes& rhs) {
        Range3f merged = lhs;
        merged.ExtendByPoint(rhs.position);
        return merged;
      });
  EXPECT_EQ(box, Range3f(kMin, kMax));
  for (auto const& vertex : box_edges) {
    EXPECT_TRUE(box.ContainsPoint(vertex.position));
  }
}

TEST(ShapeComponentUtilsTest, ShapeComponentWorks) {
  const Point3f kMin(-1.0f, 2.0f, 0.0f);
  const Point3f kMax(-0.5f, 4.0f, 3.1f);
  const Color4f kYellow(1.0f, 1.0f, 0.0f, 1.0f);

  std::unique_ptr<const ShapeComponent> box_shape =
      MakeWireframeBoxShape("wireframe_box", kMin, kMax, kYellow);
  EXPECT_EQ(box_shape->GetBoundingBox(), Range3f(kMin, kMax));
}

}  // namespace
}  // namespace component
}  // namespace seurat
