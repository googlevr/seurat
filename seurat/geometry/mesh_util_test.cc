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

#include "seurat/geometry/mesh_util.h"

#include <vector>

#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/geometry/mesh.h"
#include "seurat/geometry/triangle.h"

namespace seurat {
namespace geometry {
namespace {

using ion::math::Point3f;

TEST(MeshUtilTest, AppendTriangleFan_DegeneratePolygon) {
  Mesh mesh(1);
  std::vector<Point3f> single_position = {Point3f::Zero()};
  std::vector<Point3f> single_tex_coord = {Point3f::Zero()};

  AppendTriangleFan(single_position, {single_tex_coord}, &mesh);
  EXPECT_EQ(0, mesh.GetTriangleCount());

  std::vector<Point3f> line_segment_positions = {Point3f(1.0f, 0.0f, 0.5f),
                                                 Point3f(1.0f, 1.0f, 0.5f)};
  std::vector<Point3f> line_segment_tex_coords = {Point3f(1.0f, 0.0f, 1.0f),
                                                  Point3f(1.0f, 1.0f, 1.0f)};

  AppendTriangleFan(single_position, {single_tex_coord}, &mesh);
  EXPECT_EQ(0, mesh.GetTriangleCount());
}

TEST(MeshUtilTest, AppendTriangleFan_SingleTriangleWithTexCoords) {
  Mesh mesh(1);
  Triangle3f positions = {{
      Point3f(0.0f, 0.0f, 0.5f),  //
      Point3f(1.0f, 0.0f, 0.5f),  //
      Point3f(0.0f, 2.0f, 0.5f),
  }};

  Triangle3f tex_coords = {{
      Point3f(0.1f, 0.1f, 1.0f),  //
      Point3f(1.1f, 0.1f, 1.0f),  //
      Point3f(0.1f, 2.1f, 1.0f),
  }};

  AppendTriangleFan(positions, {tex_coords}, &mesh);
  EXPECT_EQ(1, mesh.GetTriangleCount());

  Triangle3f mesh_vertex_positions;
  Triangle3f mesh_vertex_tex_coords;

  Mesh::Triangle mesh_triangle = mesh.GetTriangles().front();
  for (int i = 0; i < 3; ++i) {
    int index = mesh_triangle[i];
    mesh_vertex_positions[i] = mesh.GetPositions()[index];
    mesh_vertex_tex_coords[i] = mesh.GetTexCoords(0)[index];
  }

  EXPECT_THAT(mesh_vertex_positions,
              ::testing::UnorderedElementsAreArray(positions));
  EXPECT_THAT(mesh_vertex_tex_coords,
              ::testing::UnorderedElementsAreArray(tex_coords));

  EXPECT_NEAR(0.0f,
              ion::math::Length(NormalFromTriangle(mesh_vertex_positions) -
                                NormalFromTriangle(positions)),
              1e-6f);
}

TEST(MeshUtilTest, AppendTriangleFan_SingleTriangleWithoutTexCoords) {
  Mesh mesh(0);
  Triangle3f positions = {{
      Point3f(0.0f, 0.0f, 0.5f),  //
      Point3f(1.0f, 0.0f, 0.5f),  //
      Point3f(0.0f, 2.0f, 0.5f),
  }};

  AppendTriangleFan(positions, {}, &mesh);
  EXPECT_EQ(1, mesh.GetTriangleCount());

  std::vector<Point3f> mesh_vertex_positions;
  std::vector<Point3f> mesh_vertex_tex_coords;
  Mesh::Triangle mesh_triangle = mesh.GetTriangles().front();
  for (const int index : mesh_triangle) {
    mesh_vertex_positions.push_back(mesh.GetPositions()[index]);
  }

  EXPECT_THAT(mesh_vertex_positions,
              ::testing::UnorderedElementsAreArray(positions));
}

TEST(MeshUtilTest, AppendTriangleFan_MultipleTriangles) {
  Mesh mesh(1);
  // Append a pentagon to the mesh:
  //
  //    *
  //   / \
  //  /   \
  // *     *
  // |     |
  // |     |
  // *-----*
  std::vector<Point3f> positions = {{
      Point3f(0.0f, 0.0f, 0.5f),  //
      Point3f(1.0f, 0.0f, 0.5f),  //
      Point3f(1.0f, 1.0f, 0.5f),  //
      Point3f(0.5f, 2.0f, 0.5f),  //
      Point3f(0.0f, 1.0f, 0.5f),  //
  }};

  // Note that 'tex_coords' are the same as 'positions' with the z-coordinate
  // dropped.
  std::vector<Point3f> tex_coords = {{
      Point3f(0.0f, 0.0f, 1.0f),  //
      Point3f(1.0f, 0.0f, 1.0f),  //
      Point3f(1.0f, 1.0f, 1.0f),  //
      Point3f(0.5f, 2.0f, 1.0f),  //
      Point3f(0.0f, 1.0f, 1.0f),  //
  }};

  AppendTriangleFan(positions, {{tex_coords}}, &mesh);
  EXPECT_EQ(3, mesh.GetTriangleCount());

  // Resulting triangles must have correct texture coordinates for each vertex.
  for (const auto& triangle : mesh.GetTriangles()) {
    for (const int index : triangle) {
      EXPECT_EQ(mesh.GetTexCoords(0)[index][0], mesh.GetPositions()[index][0]);
      EXPECT_EQ(mesh.GetTexCoords(0)[index][1], mesh.GetPositions()[index][1]);
    }
  }

  // All resulting triangles must share the first coordinate.
  for (const auto& triangle : mesh.GetTriangles()) {
    std::vector<Point3f> triangle_positions;
    std::vector<Point3f> triangle_tex_coords;
    for (const int index : triangle) {
      triangle_positions.push_back(mesh.GetPositions()[index]);
      triangle_tex_coords.push_back(mesh.GetTexCoords(0)[index]);
    }

    EXPECT_THAT(triangle_positions, ::testing::Contains(positions[0]));
    EXPECT_THAT(triangle_tex_coords, ::testing::Contains(tex_coords[0]));
  }
}

TEST(MeshUtilTest, EstimateProjectedArea_EmptyMesh) {
  Mesh mesh(0);
  EXPECT_EQ(0.0f, EstimateProjectedArea(Point3f::Zero(), mesh));
}

TEST(MeshUtilTest, EstimateProjectedArea_TessellatedPlane) {
  Mesh mesh(0);

  // Add a dense grid of triangles with z=-1 and with (x, y) in [-1.0, 1.0]^2.
  //
  // This is essentially one face of a cubemap
  const int kResolution = 50;
  for (int y = 0; y <= kResolution; ++y) {
    for (int x = 0; x <= kResolution; ++x) {
      Point3f p(x, y, 0.0f);
      p /= static_cast<float>(kResolution);
      p -= {0.5f, 0.5f, 0.0f};
      p *= 2.0f;
      p[2] = -1.0f;
      mesh.AppendVertex(p, {});
    }
  }
  for (int y = 0; y < kResolution; ++y) {
    for (int x = 0; x < kResolution; ++x) {
      int bottom_left = y * (kResolution + 1) + x;
      int bottom_right = y * (kResolution + 1) + x + 1;
      int top_left = (y + 1) * (kResolution + 1) + x;
      int top_right = (y + 1) * (kResolution + 1) + x + 1;
      mesh.AppendTriangle(
          Mesh::Triangle{{bottom_left, bottom_right, top_left}});
      mesh.AppendTriangle(Mesh::Triangle{{bottom_right, top_left, top_right}});
    }
  }

  float area_from_origin = EstimateProjectedArea(Point3f::Zero(), mesh);
  // The grid is basically one face of a cubemap, so the area should be 1/6.
  EXPECT_NEAR(1.0f / 6.0f, area_from_origin, 1e-2f);

  float area_from_near =
      EstimateProjectedArea(Point3f(0.0f, 0.0f, 100.0f - 1.0f), mesh);
  // Compute area from a point which is 2x further, so the projected area should
  // be 4x smaller.
  float area_from_far =
      EstimateProjectedArea(Point3f(0.0f, 0.0f, 200.0f - 1.0f), mesh);
  EXPECT_NEAR(4.0f, area_from_near / area_from_far, 1e-2f);
}
TEST(MeshUtilTest, ToUnindexedMesh) {
  Mesh indexed(1);
  // Points & texture coordinates for testing (power of 2's used to enable
  // float-equality comparisons without precision problems).
  const Point3f pa(2.0f, 4.0f, 8.0f);
  const Point3f pb(16.0f, 32.0f, 64.0f);
  const Point3f pc(4.0f, 8.0f, 16.0f);
  const Point3f pd(4.0f, 2.0f, 32.0f);
  const Point3f ta(2.0f, 4.0f, 1.0f);
  const Point3f tb(16.0f, 64.0f, 1.0f);
  const Point3f tc(8.0f, 16.0f, 1.0f);
  const Point3f td(8.0f, 2.0f, 1.0f);
  indexed.AppendVertex(pa, {ta});
  indexed.AppendVertex(pb, {tb});
  indexed.AppendVertex(pc, {tc});
  indexed.AppendVertex(pd, {td});

  indexed.AppendTriangle({{0, 1, 2}});
  indexed.AppendTriangle({{1, 2, 3}});
  indexed.AppendTriangle({{2, 1, 2}});

  Mesh unindexed = ToUnindexedMesh(indexed);
  for (int tri_index = 0; tri_index < unindexed.GetTriangleCount();
       ++tri_index) {
    // The unindexed triangle indices should be in order.
    Mesh::Triangle expected_triangle{
        {tri_index * 3, tri_index * 3 + 1, tri_index * 3 + 2}};
    Mesh::Triangle unindexed_triangle = unindexed.GetTriangles()[tri_index];
    EXPECT_EQ(expected_triangle, unindexed_triangle);

    // The vertex data should be the same as the corresponding triangle in the
    // indexed mesh.
    Mesh::Triangle indexed_triangle = indexed.GetTriangles()[tri_index];
    for (int i = 0; i < 3; ++i) {
      EXPECT_EQ(indexed.GetPositions()[indexed_triangle[i]],
                unindexed.GetPositions()[unindexed_triangle[i]]);
      EXPECT_EQ(indexed.GetTexCoords(0)[indexed_triangle[i]],
                unindexed.GetTexCoords(0)[unindexed_triangle[i]]);
    }
  }
}

TEST(MeshUtilTest, SplitMeshByVertexCount) {
  const Point3f pa(2.0f, 4.0f, 8.0f);
  const Point3f pb(16.0f, 32.0f, 64.0f);
  const Point3f pc(4.0f, 8.0f, 16.0f);
  const Point3f pd(4.0f, 2.0f, 32.0f);
  const Point3f ta(2.0f, 4.0f, 1.0f);
  const Point3f tb(16.0f, 64.0f, 1.0f);
  const Point3f tc(8.0f, 16.0f, 1.0f);
  const Point3f td(8.0f, 2.0f, 1.0f);
  Mesh original(1);
  original.AppendVertex(pa, {{ta}});
  original.AppendVertex(pb, {{tb}});
  original.AppendVertex(pc, {{tc}});
  original.AppendVertex(pd, {{td}});

  // Add triangles such that we expect the mesh to be split as follows:
  // Mesh 0, contains vertices {0, 1, 2}.
  original.AppendTriangle({{0, 1, 2}});

  // Mesh 1, contains vertices {1, 2, 3}
  original.AppendTriangle({{1, 2, 3}});
  original.AppendTriangle({{2, 1, 2}});
  original.AppendTriangle({{2, 2, 2}});

  // Mesh 2, contains vertices {0, 1, 3}
  original.AppendTriangle({{0, 0, 0}});
  original.AppendTriangle({{3, 3, 3}});
  original.AppendTriangle({{3, 3, 3}});
  original.AppendTriangle({{1, 1, 3}});
  original.AppendTriangle({{1, 1, 3}});

  // Mesh 3, contains vertices {0, 2}
  original.AppendTriangle({{0, 2, 2}});

  // Mesh 4, contains vertices {1, 2, 3}
  original.AppendTriangle({{1, 2, 3}});

  std::vector<Mesh> noop_split =
      SplitMeshByVertexCount(original, original.GetVertexCount());
  EXPECT_EQ(1, noop_split.size());
  EXPECT_EQ(original.GetPositions(), noop_split.front().GetPositions());
  EXPECT_EQ(original.GetTexCoords(0), noop_split.front().GetTexCoords(0));
  EXPECT_EQ(original.GetTriangles(), noop_split.front().GetTriangles());

  std::vector<Mesh> split = SplitMeshByVertexCount(original, 3);
  EXPECT_EQ(5, split.size());

  EXPECT_EQ(1, split[0].GetTextureCount());
  EXPECT_EQ(1, split[0].GetTriangleCount());
  EXPECT_THAT(split[0].GetPositions(),
              ::testing::UnorderedElementsAreArray(
                  std::vector<Point3f>{original.GetPositions()[0],  //
                                       original.GetPositions()[1],  //
                                       original.GetPositions()[2]}));
  EXPECT_THAT(split[0].GetTexCoords(0),
              ::testing::UnorderedElementsAreArray(
                  std::vector<Point3f>{original.GetTexCoords(0)[0],  //
                                       original.GetTexCoords(0)[1],  //
                                       original.GetTexCoords(0)[2]}));

  EXPECT_EQ(1, split[1].GetTextureCount());
  EXPECT_EQ(3, split[1].GetTriangleCount());
  EXPECT_THAT(split[1].GetPositions(),
              ::testing::UnorderedElementsAreArray(
                  std::vector<Point3f>{original.GetPositions()[1],  //
                                       original.GetPositions()[2],  //
                                       original.GetPositions()[3]}));
  EXPECT_THAT(split[1].GetTexCoords(0),
              ::testing::UnorderedElementsAreArray(
                  std::vector<Point3f>{original.GetTexCoords(0)[1],  //
                                       original.GetTexCoords(0)[2],  //
                                       original.GetTexCoords(0)[3]}));

  EXPECT_EQ(1, split[2].GetTextureCount());
  EXPECT_EQ(5, split[2].GetTriangleCount());
  EXPECT_THAT(split[2].GetPositions(),
              ::testing::UnorderedElementsAreArray(
                  std::vector<Point3f>{original.GetPositions()[0],  //
                                       original.GetPositions()[1],  //
                                       original.GetPositions()[3]}));
  EXPECT_THAT(split[2].GetTexCoords(0),
              ::testing::UnorderedElementsAreArray(
                  std::vector<Point3f>{original.GetTexCoords(0)[0],  //
                                       original.GetTexCoords(0)[1],  //
                                       original.GetTexCoords(0)[3]}));

  EXPECT_EQ(1, split[3].GetTextureCount());
  EXPECT_EQ(1, split[3].GetTriangleCount());
  EXPECT_THAT(split[3].GetPositions(),
              ::testing::UnorderedElementsAreArray(
                  std::vector<Point3f>{original.GetPositions()[0],  //
                                       original.GetPositions()[2]}));
  EXPECT_THAT(split[3].GetTexCoords(0),
              ::testing::UnorderedElementsAreArray(
                  std::vector<Point3f>{original.GetTexCoords(0)[0],  //
                                       original.GetTexCoords(0)[2]}));

  EXPECT_EQ(1, split[4].GetTextureCount());
  EXPECT_EQ(1, split[4].GetTriangleCount());
  EXPECT_THAT(split[4].GetPositions(),
              ::testing::UnorderedElementsAreArray(
                  std::vector<Point3f>{original.GetPositions()[1],  //
                                       original.GetPositions()[2],  //
                                       original.GetPositions()[3]}));
  EXPECT_THAT(split[4].GetTexCoords(0),
              ::testing::UnorderedElementsAreArray(
                  std::vector<Point3f>{original.GetTexCoords(0)[1],  //
                                       original.GetTexCoords(0)[2],  //
                                       original.GetTexCoords(0)[3]}));
}

}  // namespace
}  // namespace geometry
}  // namespace seurat
