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
#include <limits>

#include "ion/math/matrixutils.h"
#include "ion/math/vectorutils.h"
#include "seurat/base/color.h"
#include "seurat/geometry/mesh.h"
#include "seurat/image/image.h"
#include "seurat/image/image_util.h"

using ion::math::Matrix3f;
using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Vector2f;
using ion::math::Vector3f;
using seurat::base::Color4ui;
using seurat::geometry::Mesh;
using seurat::image::Image4ui8;

namespace seurat {
namespace mesh {

// static
std::unique_ptr<const MeshComponent> MeshComponentUtil::CreateCube(
    ion::gfx::ImagePtr cube_face_texture) {
  std::vector<MeshComponent::VertexAttributes> attribute_buf;
  std::vector<uint32> index_buf;

  std::array<Point3f, 4> base_face_verts = {{{-0.5f, -0.5f, -0.5f},
                                             {0.5f, -0.5f, -0.5f},
                                             {-0.5f, 0.5f, -0.5f},
                                             {0.5f, 0.5f, -0.5f}}};
  std::array<Point3f, 4> tex_coords = {{{0.0f, 0.0f, 1.0f},
                                        {1.0f, 0.0f, 1.0f},
                                        {0.0f, 1.0f, 1.0f},
                                        {1.0f, 1.0f, 1.0f}}};

  for (int i = 0; i < 6; ++i) {
    // Create a rotation matrix to place base_face_verts into a new "face space"
    // which is rotated to position a face along the xy-plane onto the plane
    // consisting of face_x and face_y, with normal face_z.
    Vector3f face_x = Vector3f::Zero();
    face_x[i % 3] = 1.0f;
    Vector3f face_y = Vector3f::Zero();
    face_y[(i + 1) % 3] = 1.0f;
    if (i >= 3) {
      std::swap(face_x, face_y);
    }
    Vector3f face_z = ion::math::Cross(face_x, face_y);
    Matrix3f rotation(face_x[0], face_y[0], face_z[0],  //
                      face_x[1], face_y[1], face_z[1],  //
                      face_x[2], face_y[2], face_z[2]);
    // Ensure that we are not reflecting the base face.
    CHECK_GT(ion::math::Determinant(rotation), 0.0f);

    int base_index = attribute_buf.size();
    for (int j = 0; j < 4; ++j) {
      const Point3f& vert = base_face_verts[j];
      const auto& tex_coord = tex_coords[j];
      attribute_buf.push_back({rotation * vert, tex_coord});
    }
    index_buf.push_back(base_index + 0);
    index_buf.push_back(base_index + 1);
    index_buf.push_back(base_index + 3);

    index_buf.push_back(base_index + 0);
    index_buf.push_back(base_index + 3);
    index_buf.push_back(base_index + 2);
  }

  return MeshComponent::Create("cube", cube_face_texture,
                               std::move(attribute_buf), std::move(index_buf));
}

Mesh MeshComponentUtil::ToMesh(const MeshComponent& component) {
  Mesh mesh(1);
  for (const auto& vertex : component.GetAttributeBuffer()) {
    mesh.AppendVertex(vertex.position, {vertex.tex_coord});
  }
  CHECK_LE(component.GetAttributeBuffer().size(),
           static_cast<size_t>(std::numeric_limits<int>::max()));
  const int num_triangles = component.GetIndexBuffer().size() / 3;
  for (int t = 0; t < num_triangles; ++t) {
    Mesh::Triangle triangle;
    triangle[0] = component.GetIndexBuffer()[3 * t];
    triangle[1] = component.GetIndexBuffer()[3 * t + 1];
    triangle[2] = component.GetIndexBuffer()[3 * t + 2];
    mesh.AppendTriangle(triangle);
  }
  return mesh;
}

std::unique_ptr<const MeshComponent> MeshComponentUtil::FromMesh(
    const Mesh& mesh, ion::gfx::ImagePtr texture_atlas) {
  CHECK_EQ(1, mesh.GetTextureCount());
  std::vector<MeshComponent::VertexAttributes> attributes;
  attributes.reserve(mesh.GetVertexCount());
  for (int i = 0; i < mesh.GetVertexCount(); ++i) {
    attributes.push_back({mesh.GetPositions()[i], mesh.GetTexCoords(0)[i]});
  }

  std::vector<uint32> index_buf;
  index_buf.reserve(mesh.GetTriangleCount() * 3);
  for (const Mesh::Triangle& triangle : mesh.GetTriangles()) {
    for (int i = 0; i < 3; ++i) {
      DCHECK_GE(triangle[i], 0);
      DCHECK_LT(triangle[i], mesh.GetVertexCount());
      index_buf.push_back(static_cast<uint32>(triangle[i]));
    }
  }

  return MeshComponent::Create("mesh", texture_atlas, std::move(attributes),
                               std::move(index_buf));
}

std::unique_ptr<const MeshComponent>
MeshComponentUtil::SortTrianglesByDistanceToEye(const MeshComponent& original,
                                                const ion::math::Point3f& eye,
                                                SortingOrder sorting_order) {
  struct Triangle {
    bool operator<(const Triangle& rhs) const {
      return distance_squared < rhs.distance_squared;
    }

    float distance_squared;
    std::array<uint32, 3> indices;
  };

  const std::vector<MeshComponent::VertexAttributes>& attrs =
      original.GetAttributeBuffer();
  const std::vector<uint32>& index_buf = original.GetIndexBuffer();

  std::vector<Triangle> triangles;
  CHECK_EQ(index_buf.size() % 3, 0);
  for (int i = 0; i < index_buf.size(); i += 3) {
    Point3f center = (attrs[index_buf[i + 0]].position +  //
                      attrs[index_buf[i + 1]].position +  //
                      attrs[index_buf[i + 2]].position) /
                     3.0f;
    triangles.push_back({
        ion::math::LengthSquared(center - eye),
        {{
            index_buf[i + 0],  //
            index_buf[i + 1],  //
            index_buf[i + 2]   //
        }},
    });
  }

  std::vector<uint32> final_index_buf;
  final_index_buf.reserve(index_buf.size());
  std::sort(triangles.begin(), triangles.end());
  switch (sorting_order) {
    case SortingOrder::kFrontToBack:
      for (const Triangle& tri : triangles) {
        final_index_buf.push_back(tri.indices[0]);
        final_index_buf.push_back(tri.indices[1]);
        final_index_buf.push_back(tri.indices[2]);
      }
      break;
    case SortingOrder::kBackToFront:
      for (auto tri = triangles.rbegin(); tri != triangles.rend(); ++tri) {
        final_index_buf.push_back(tri->indices[0]);
        final_index_buf.push_back(tri->indices[1]);
        final_index_buf.push_back(tri->indices[2]);
      }
      break;
  }
  return MeshComponent::Create(original.GetLabel(), original.GetTextureAtlas(),
                               attrs, std::move(final_index_buf));
}

int MeshComponentUtil::CountTransparentTexelEdges(const MeshComponent& mesh) {
  const ion::gfx::ImagePtr& texture = mesh.GetTextureAtlas();
  Image4ui8 image = image::ConvertIonImageToSeuratImage<Image4ui8>(texture);
  int edges = 0;
  for (int y = 0; y < image.Height(); ++y) {
    for (int x = 0; x < image.Width(); ++x) {
      bool cur_pixel_opaque = image.At(x, y)[3] > 0;
      if (y + 1 < texture->GetHeight()) {
        bool top_pixel_opaque = image.At(x, y + 1)[3] > 0;
        // bitwise-xor works because false and true are guaranteed to be cast to
        // 0 and 1, respectively.
        edges += cur_pixel_opaque ^ top_pixel_opaque;
      }
      if (x + 1 < texture->GetWidth()) {
        bool right_pixel_opaque = image.At(x + 1, y)[3] > 0;
        edges += cur_pixel_opaque ^ right_pixel_opaque;
      }
    }
  }
  return edges;
}

}  // namespace mesh
}  // namespace seurat
