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

#include "seurat/geometry/mesh.h"

#include <vector>

#include "ion/math/matrixutils.h"
#include "ion/math/transformutils.h"
#include "gtest/gtest.h"

namespace seurat {
namespace geometry {

using ion::math::Point3f;
using Triangle = Mesh::Triangle;

TEST(Mesh, Construction) {
  Mesh mesh(1);

  EXPECT_EQ(1, mesh.GetTextureCount());
  EXPECT_EQ(0, mesh.GetVertexCount());
  EXPECT_EQ(0, mesh.GetTriangleCount());

  const std::vector<Point3f> positions{
      Point3f(0.0f, 0.0f, 0.0f), Point3f(1.0f, 0.0f, 0.0f),
      Point3f(1.0f, 1.0f, 0.0f), Point3f(0.0f, 1.0f, 0.0f)};
  const std::vector<std::vector<Point3f>> tex_coords{
      {Point3f(0.0f, 0.0f, 1.0f), Point3f(1.0f, 0.0f, 1.0f),
       Point3f(1.0f, 1.0f, 1.0f), Point3f(0.0f, 1.0f, 1.0f)}};
  std::vector<Triangle> triangles{{{0, 1, 2}}, {{2, 3, 0}}};

  for (int i = 0; i < 4; ++i) {
    mesh.AppendVertex(positions[i], {{tex_coords[0][i]}});
  }
  EXPECT_EQ(4, mesh.GetVertexCount());
  EXPECT_EQ(absl::Span<const Point3f>(positions), mesh.GetPositions());
  for (int i = 0; i < mesh.GetTextureCount(); ++i) {
    EXPECT_EQ(absl::Span<const Point3f>(tex_coords[i]), mesh.GetTexCoords(i));
  }

  for (const Triangle& triangle : triangles) {
    mesh.AppendTriangle(triangle);
  }
  EXPECT_EQ(2, mesh.GetTriangleCount());
  EXPECT_EQ(absl::Span<const Triangle>(triangles), mesh.GetTriangles());
}

TEST(Mesh, AppendMesh) {
  Mesh mesh_a(1);

  EXPECT_EQ(1, mesh_a.GetTextureCount());

  std::vector<Point3f> expected_positions;
  std::vector<std::vector<Point3f>> expected_tex_coords;
  expected_tex_coords.resize(mesh_a.GetTextureCount());
  std::vector<Triangle> expected_triangles;

  const std::vector<Point3f> positions_a{Point3f(0.0f, 0.0f, 0.0f),
                                         Point3f(1.0f, 0.0f, 0.0f),
                                         Point3f(1.0f, 1.0f, 0.0f)};
  const std::vector<std::vector<Point3f>> tex_coords_a{
      {Point3f(0.0f, 0.0f, 1.0f), Point3f(1.0f, 0.0f, 1.0f),
       Point3f(1.0f, 1.0f, 1.0f)}};
  EXPECT_EQ(3, positions_a.size());
  EXPECT_EQ(mesh_a.GetTextureCount(), tex_coords_a.size());

  for (int i = 0; i < 3; ++i) {
    mesh_a.AppendVertex(positions_a[i], {{tex_coords_a[0][i]}});
  }
  expected_positions = positions_a;
  expected_tex_coords = tex_coords_a;
  std::vector<Triangle> triangles_a{Triangle{{0, 1, 2}}};
  for (const Triangle& triangle : triangles_a) {
    mesh_a.AppendTriangle(triangle);
    expected_triangles.push_back(triangle);
  }

  EXPECT_EQ(3, mesh_a.GetVertexCount());
  EXPECT_EQ(1, mesh_a.GetTriangleCount());

  Mesh mesh_b(1);
  const std::vector<Point3f> positions_b{Point3f(0.0f, 0.0f, 1.0f),
                                         Point3f(1.0f, 0.0f, 1.0f),
                                         Point3f(1.0f, 1.0f, 1.0f)};
  const std::vector<std::vector<Point3f>> tex_coords_b{
      {Point3f(0.0f, 0.0f, 1.0f), Point3f(1.0f, 0.0f, 1.0f),
       Point3f(1.0f, 1.0f, 1.0f)}};

  for (int i = 0; i < 3; ++i) {
    expected_positions.push_back(positions_b[i]);
    mesh_b.AppendVertex(positions_b[i], {{tex_coords_b[0][i]}});
    expected_tex_coords[0].push_back(tex_coords_b[0][i]);
  }
  std::vector<Triangle> triangles_b{Triangle{{0, 1, 2}}};
  for (const Triangle& triangle : triangles_b) {
    mesh_b.AppendTriangle(triangle);
    expected_triangles.push_back(
        Triangle{{triangle[0] + mesh_a.GetVertexCount(),
                  triangle[1] + mesh_a.GetVertexCount(),
                  triangle[2] + mesh_a.GetVertexCount()}});
  }

  EXPECT_EQ(1, mesh_b.GetTextureCount());
  EXPECT_EQ(3, mesh_b.GetVertexCount());
  EXPECT_EQ(1, mesh_b.GetTriangleCount());

  // Append to empty mesh.
  Mesh mesh(1);
  mesh.AppendMesh(mesh_a);
  mesh.AppendMesh(mesh_b);
  EXPECT_EQ(1, mesh.GetTextureCount());
  EXPECT_EQ(6, mesh.GetVertexCount());
  EXPECT_EQ(2, mesh.GetTriangleCount());
  EXPECT_EQ(absl::Span<const Point3f>(expected_positions), mesh.GetPositions());
  for (int i = 0; i < mesh.GetTextureCount(); ++i) {
    EXPECT_EQ(absl::Span<const Point3f>(expected_tex_coords[i]),
              mesh.GetTexCoords(i));
  }
  for (int i = 0; i < mesh.GetTriangleCount(); ++i) {
    EXPECT_EQ(expected_triangles[i], mesh.GetTriangles()[i]);
  }

  // Append to non-empty mesh.
  mesh_a.AppendMesh(mesh_b);
  EXPECT_EQ(6, mesh_a.GetVertexCount());
  EXPECT_EQ(2, mesh_a.GetTriangleCount());
  EXPECT_EQ(absl::Span<const Point3f>(expected_positions),
            mesh_a.GetPositions());
  for (int i = 0; i < mesh_a.GetTextureCount(); ++i) {
    EXPECT_EQ(absl::Span<const Point3f>(expected_tex_coords[i]),
              mesh_a.GetTexCoords(i));
  }
  for (int i = 0; i < mesh_a.GetTriangleCount(); ++i) {
    EXPECT_EQ(expected_triangles[i], mesh_a.GetTriangles()[i]);
  }
}

TEST(Mesh, Transform) {
  const float kLeft = -10.0f;
  const float kRight = 10.0f;
  const float kBottom = -10.0f;
  const float kTop = 10.0f;
  const float kNear = 10.0f;
  const float kFar = 100.0f;
  const ion::math::Matrix4f perspective_matrix =
      ion::math::PerspectiveMatrixFromFrustum(kLeft, kRight, kBottom, kTop,
                                              kNear, kFar);

  const ion::math::Matrix4f inverse_perspective_matrix =
      ion::math::Inverse(perspective_matrix);

  Mesh original_mesh(1);
  original_mesh.AppendVertex(Point3f(0.0f, 0.0f, -10.0f),
                             {{Point3f(0.0f, 0.0f, 1.0f)}});
  original_mesh.AppendVertex(Point3f(1.0f, 0.0f, -10.0f),
                             {{Point3f(1.0f, 0.0f, 1.0f)}});
  original_mesh.AppendVertex(Point3f(1.0f, 1.0f, -10.0f),
                             {{Point3f(1.0f, 1.0f, 1.0f)}});
  original_mesh.AppendTriangle({{0, 1, 2}});

  Mesh transformed = original_mesh;
  transformed.TransformPositions(perspective_matrix);

  EXPECT_NE(original_mesh.GetPositions(), transformed.GetPositions());
  EXPECT_EQ(
      ion::math::ProjectPoint(perspective_matrix, Point3f(1.0f, 0.0f, -10.0f)),
      transformed.GetPositions()[1]);
  EXPECT_EQ(original_mesh.GetTexCoords(0), transformed.GetTexCoords(0));
  EXPECT_EQ(original_mesh.GetTriangles(), transformed.GetTriangles());

  Mesh transformed2 = transformed;
  transformed2.TransformPositions(inverse_perspective_matrix);
  EXPECT_EQ(original_mesh.GetPositions(), transformed2.GetPositions());
  EXPECT_EQ(original_mesh.GetTexCoords(0), transformed.GetTexCoords(0));
  EXPECT_EQ(original_mesh.GetTriangles(), transformed2.GetTriangles());
}

}  // namespace geometry
}  // namespace seurat
