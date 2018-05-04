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

#include "seurat/mesh/mesh_component_util.h"

#include <array>
#include <random>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/geometry/mesh.h"

using ion::math::Point3f;
using ion::math::Vector3f;
using seurat::geometry::Mesh;

namespace seurat {
namespace mesh {
namespace {

TEST(MeshComponentUtilTest, TestCreateCube) {
  ion::gfx::ImagePtr texture =
      base::CreateImage(ion::gfx::Image::kRgba8888, {16, 16});
  const auto cube = MeshComponentUtil::CreateCube(texture);
  const auto attrs = cube->GetAttributeBuffer();
  const auto tri_indices = cube->GetIndexBuffer();

  // A cube should have 6 faces, each with 2 triangles consisting of 3 vertices.
  EXPECT_EQ(6 * 2 * 3, tri_indices.size());

  // The cube should be unit-sized and origin-centered.
  EXPECT_EQ(Vector3f(1.0f, 1.0f, 1.0f), cube->GetBoundingBox().GetSize());
  EXPECT_EQ(Point3f::Zero(), cube->GetBoundingBox().GetCenter());
}

TEST(MeshComponentUtilTest, TestCubeFaceOrientation) {
  ion::gfx::ImagePtr texture =
      base::CreateImage(ion::gfx::Image::kRgba8888, {16, 16});
  const auto cube = MeshComponentUtil::CreateCube(texture);
  const auto attrs = cube->GetAttributeBuffer();
  const auto tri_indices = cube->GetIndexBuffer();

  for (int tri_idx = 0; tri_idx < tri_indices.size() / 3; ++tri_idx) {
    const auto& attr0 = attrs[tri_indices[tri_idx * 3 + 0]];
    const auto& attr1 = attrs[tri_indices[tri_idx * 3 + 1]];
    const auto& attr2 = attrs[tri_indices[tri_idx * 3 + 2]];

    const Point3f center =
        (attr0.position + attr1.position + attr2.position) / 3.0f;
    // If 0->1->2 form a counter-clockwise triangle, then this normal should
    // face the outside of the cube, in the opposite direction of the
    // triangle-center->origin vector.
    const Vector3f normal = ion::math::Cross(attr2.position - attr1.position,
                                             attr0.position - attr1.position);

    EXPECT_GT(0.0f, ion::math::Dot(normal, center - Point3f::Zero()));
  }
}

TEST(MeshComponentUtilTest, ToMesh) {
  ion::gfx::ImagePtr texture =
      base::CreateImage(ion::gfx::Image::kRgba8888, {2, 2});
  // Generate a 1x1 square mesh in the xy-plane.
  const auto component = MeshComponent::Create(
      "1x1_square", texture,
      {
          {Point3f(0.0f, 0.0f, 0.0f), Point3f(0.0f, 0.0f, 1.0f)},
          {Point3f(1.0f, 0.0f, 0.0f), Point3f(1.0f, 0.0f, 1.0f)},
          {Point3f(0.0f, 1.0f, 0.0f), Point3f(0.0f, 1.0f, 1.0f)},
          {Point3f(1.0f, 1.0f, 0.0f), Point3f(1.0f, 1.0f, 1.0f)},
      },
      {0, 1, 3, 0, 3, 2});

  Mesh mesh = MeshComponentUtil::ToMesh(*component);

  EXPECT_EQ(4, mesh.GetVertexCount());
  EXPECT_EQ(2, mesh.GetTriangleCount());
  EXPECT_EQ((std::vector<Mesh::Triangle>{Mesh::Triangle{{0, 1, 3}},
                                         Mesh::Triangle{{0, 3, 2}}}),
            mesh.GetTriangles());
  EXPECT_EQ((std::vector<Point3f>{
                Point3f(0.0f, 0.0f, 0.0f), Point3f(1.0f, 0.0f, 0.0f),
                Point3f(0.0f, 1.0f, 0.0f), Point3f(1.0f, 1.0f, 0.0f)}),
            mesh.GetPositions());
  EXPECT_EQ((std::vector<Point3f>{
                Point3f(0.0f, 0.0f, 1.0f), Point3f(1.0f, 0.0f, 1.0f),
                Point3f(0.0f, 1.0f, 1.0f), Point3f(1.0f, 1.0f, 1.0f)}),
            mesh.GetTexCoords(0));
}

TEST(MeshComponentUtilTest, FromMesh) {
  ion::gfx::ImagePtr texture =
      base::CreateImage(ion::gfx::Image::kRgba8888, {2, 2});
  Mesh mesh(1);
  mesh.AppendVertex(Point3f(0.0f, 0.0f, 0.0f), {Point3f(0.0f, 0.0f, 1.0f)});
  mesh.AppendVertex(Point3f(1.0f, 0.0f, 0.0f), {Point3f(1.0f, 0.0f, 1.0f)});
  mesh.AppendVertex(Point3f(0.0f, 1.0f, 0.0f), {Point3f(0.0f, 1.0f, 1.0f)});
  mesh.AppendVertex(Point3f(1.0f, 1.0f, 0.0f), {Point3f(1.0f, 1.0f, 1.0f)});
  mesh.AppendTriangle(Mesh::Triangle{{0, 1, 3}});
  mesh.AppendTriangle(Mesh::Triangle{{0, 3, 2}});

  const auto component = MeshComponentUtil::FromMesh(mesh, texture);
  EXPECT_EQ(4U, component->GetAttributeBuffer().size());
  EXPECT_EQ(6U, component->GetIndexBuffer().size());
  std::vector<MeshComponent::VertexAttributes> expected_attribute_buffer{
      {Point3f(0.0f, 0.0f, 0.0f), Point3f(0.0f, 0.0f, 1.0f)},
      {Point3f(1.0f, 0.0f, 0.0f), Point3f(1.0f, 0.0f, 1.0f)},
      {Point3f(0.0f, 1.0f, 0.0f), Point3f(0.0f, 1.0f, 1.0f)},
      {Point3f(1.0f, 1.0f, 0.0f), Point3f(1.0f, 1.0f, 1.0f)}};
  EXPECT_EQ(expected_attribute_buffer, component->GetAttributeBuffer());
  std::vector<uint32> expected_index_buffer{{0, 1, 3, 0, 3, 2}};
  EXPECT_EQ(expected_index_buffer, component->GetIndexBuffer());
}

TEST(MeshComponentUtilTest, SortTrianglesByDistanceToEye) {
  ion::gfx::ImagePtr texture =
      base::CreateImage(ion::gfx::Image::kRgba8888, {16, 16});
  // Two parallel triangles on the z=1 and z=2 planes.
  const auto original_mesh = MeshComponent::Create(
      "parallel_triangles", texture,
      {
          {Point3f(0.0f, 0.0f, 1.0f), Point3f(0.0f, 0.0f, 1.0f)},
          {Point3f(1.0f, 0.0f, 1.0f), Point3f(1.0f, 0.0f, 1.0f)},
          {Point3f(0.0f, 1.0f, 1.0f), Point3f(0.0f, 1.0f, 1.0f)},
          {Point3f(0.0f, 0.0f, 2.0f), Point3f(0.0f, 0.0f, 1.0f)},
          {Point3f(1.0f, 0.0f, 2.0f), Point3f(1.0f, 0.0f, 1.0f)},
          {Point3f(0.0f, 1.0f, 2.0f), Point3f(0.0f, 1.0f, 1.0f)},
      },
      {0, 1, 2, 3, 4, 5});

  const auto sorted_increasing =
      MeshComponentUtil::SortTrianglesByDistanceToEye(
          *original_mesh, Point3f(0.0f, 0.0f, 0.0f),
          MeshComponentUtil::SortingOrder::kFrontToBack);
  {
    const std::vector<uint32>& index_buf = sorted_increasing->GetIndexBuffer();
    const std::vector<MeshComponent::VertexAttributes>& attrs =
        sorted_increasing->GetAttributeBuffer();
    // Verify that the first triangle has the correct z-values for its
    // vertices.
    EXPECT_EQ(1.0f, attrs[index_buf[0]].position[2]);
    EXPECT_EQ(1.0f, attrs[index_buf[1]].position[2]);
    EXPECT_EQ(1.0f, attrs[index_buf[2]].position[2]);
    // Verify that the second triangle has the correct z-values for its
    // vertices.
    EXPECT_EQ(2.0f, attrs[index_buf[3]].position[2]);
    EXPECT_EQ(2.0f, attrs[index_buf[4]].position[2]);
    EXPECT_EQ(2.0f, attrs[index_buf[5]].position[2]);
  }

  const auto sorted_decreasing =
      MeshComponentUtil::SortTrianglesByDistanceToEye(
          *original_mesh, Point3f(0.0f, 0.0f, 0.0f),
          MeshComponentUtil::SortingOrder::kBackToFront);
  {
    const std::vector<uint32>& index_buf = sorted_decreasing->GetIndexBuffer();
    const std::vector<MeshComponent::VertexAttributes>& attrs =
        sorted_decreasing->GetAttributeBuffer();
    // Verify that the first triangle has the correct z-values for its
    // vertices.
    EXPECT_EQ(2.0f, attrs[index_buf[0]].position[2]);
    EXPECT_EQ(2.0f, attrs[index_buf[1]].position[2]);
    EXPECT_EQ(2.0f, attrs[index_buf[2]].position[2]);
    // Verify that the second triangle has the correct z-values for its
    // vertices.
    EXPECT_EQ(1.0f, attrs[index_buf[3]].position[2]);
    EXPECT_EQ(1.0f, attrs[index_buf[4]].position[2]);
    EXPECT_EQ(1.0f, attrs[index_buf[5]].position[2]);
  }
}

TEST(MeshComponentUtilTest, CountTransparentTexelEdges) {
  std::array<std::array<char, 5>, 4> test_alpha_matte = {{
      {{0, 0, 1, 0, 0}},  //
      {{0, 1, 1, 1, 0}},  //
      {{0, 1, 1, 0, 0}},  //
      {{0, 1, 1, 0, 1}},  //
  }};
  int expected_result = 13;

  // Generate a texture based on our test_alpha_matte.
  ion::gfx::ImagePtr texture =
      base::CreateImage(ion::gfx::Image::kRgba8888, {5, 4});
  const auto data_container = texture->GetData();
  uint8* texture_data = data_container->GetMutableData<uint8>();
  std::mt19937 prng;
  std::uniform_int_distribution<int> opaque_alpha_dist(1, 255);
  for (int y = 0; y < texture->GetHeight(); ++y) {
    for (int x = 0; x < texture->GetWidth(); ++x) {
      uint8 alpha;
      if (test_alpha_matte[y][x] == 1) {
        alpha = opaque_alpha_dist(prng);
      } else {
        alpha = 0;
      }
      int pixel_index = y * texture->GetWidth() + x;
      texture_data[pixel_index * 4 + 3] = alpha;
    }
  }

  const auto mesh = MeshComponent::Create("mesh", texture, {}, {});
  int result = MeshComponentUtil::CountTransparentTexelEdges(*mesh);
  EXPECT_EQ(expected_result, result);
}

}  // namespace
}  // namespace mesh
}  // namespace seurat
