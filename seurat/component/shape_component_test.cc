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

#include <memory>
#include <vector>

#include "gtest/gtest.h"
#include "seurat/base/structured_io.h"
#include "seurat/component/component_tests_utils.h"
#include "seurat/component/group_component.h"

namespace seurat {
namespace component {
namespace {

TEST(ShapeComponentTest, Construction) {
  const base::Color4f kRed{1.0f, 0.0f, 0.0f, 1.0};
  std::vector<ShapeComponent::VertexAttributes> quad_vertices = {
      {{0.0f, 0.0f, -5.0f}, kRed},    //
      {{10.0f, 0.0f, -5.0f}, kRed},   //
      {{10.0f, 10.0f, -5.0f}, kRed},  //
      {{0.0f, 10.0f, -5.0f}, kRed}    //
  };
  std::vector<uint32> quad_indices{0, 1, 2, 0, 2, 3};
  ShapeComponent triangle_shape("triangulated_quad",
                                ion::gfx::Shape::kTriangles, quad_vertices,
                                quad_indices);
  EXPECT_EQ(triangle_shape.GetLabel(), "triangulated_quad");
  EXPECT_GT(triangle_shape.GetVertices().size(), 0);
  EXPECT_GT(triangle_shape.GetIndices().size(), 0);
  EXPECT_EQ(triangle_shape.GetPrimitiveType(), ion::gfx::Shape::kTriangles);
}

TEST(ShapeComponentTest, TestIo) {
  std::vector<ShapeComponent::VertexAttributes> line_vertices = {
      {{-10.0f, 0.0f, -5.0f}, {1.0f, 0.0f, 0.0f, 1.0f}},  //
      {{10.0f, 2.0f, -5.0f}, {0.0f, 0.4f, 0.0f, 1.0f}}    //
  };
  ShapeComponent line_shape("one_line", ion::gfx::Shape::kLines, line_vertices);

  std::vector<ShapeComponent::VertexAttributes> triangle_vertices = {
      {{0.0f, 0.0f, -5.0f}, {0.0f, 0.0f, 0.0f, 1.0f}},   //
      {{10.0f, 0.0f, -5.0f}, {0.5f, 0.5f, 5.0f, 1.0f}},  //
      {{0.0f, 10.0f, -5.0f}, {1.0f, 1.0f, 1.0f, 0.5f}}   //
  };
  ShapeComponent triangle_shape("one_triangle", ion::gfx::Shape::kTriangles,
                                triangle_vertices);

  EXPECT_NE(line_shape, triangle_shape);

  ComponentCommonTests::TestIo(line_shape, triangle_shape);
}

TEST(ShapeComponentTest, TestWithIndices) {
  std::vector<ShapeComponent::VertexAttributes> line_vertices = {
      {{-10.0f, 0.0f, -5.0f}, {1.0f, 0.0f, 0.0f, 1.0f}},  //
      {{10.0f, 2.0f, -5.0f}, {0.0f, 0.25f, 1.0f, 1.0f}}   //
  };
  std::vector<uint32> line_indices{2, 1, 0};
  ShapeComponent line_shape("one_line", ion::gfx::Shape::kLines,
                            std::move(line_vertices), std::move(line_indices));

  const base::Color4f kRed{1.0f, 0.0f, 0.0f, 1.0};
  std::vector<ShapeComponent::VertexAttributes> quad_vertices = {
      {{0.0f, 0.0f, -5.0f}, kRed},    //
      {{10.0f, 0.0f, -5.0f}, kRed},   //
      {{10.0f, 10.0f, -5.0f}, kRed},  //
      {{0.0f, 10.0f, -5.0f}, kRed}    //
  };
  std::vector<uint32> quad_indices{0, 1, 2, 0, 2, 3};
  ShapeComponent triangle_shape("triangulated_quad",
                                ion::gfx::Shape::kTriangles, quad_vertices,
                                quad_indices);
  EXPECT_NE(line_shape, triangle_shape);

  ComponentCommonTests::TestIo(line_shape, triangle_shape);
}

}  // namespace
}  // namespace component
}  // namespace seurat
