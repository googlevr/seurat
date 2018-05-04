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

#ifndef VR_SEURAT_GEOMETRY_MESH_H_
#define VR_SEURAT_GEOMETRY_MESH_H_

#include <array>
#include <vector>

#include "ion/base/logging.h"
#include "ion/math/matrix.h"
#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/base/util.h"

namespace seurat {
namespace geometry {

// A simple indexed triangle set that allows for zero or more sets of texture
// coordinates. Vertex positions, vertex texture coordinates and triangles are
// stored contiguously in std::vectors. Each triangle is defined by three
// indices into the positions/texture coordinates vector. Each vertex has a
// position and a texture coordinate. The vertex attributes cannot be indexed
// separately (e.g. sharing the position, but having different texture
// coordinates is not possible).
class Mesh {
 public:
  using Triangle = std::array<int, 3>;

  // Default constructor, no texture coordinates.
  Mesh() : Mesh(0) {}

  explicit Mesh(int texture_count) {
    CHECK_LE(0, texture_count);
    tex_coords_.resize(texture_count);
  }

  // Returns the number of texture coordinate sets.
  int GetTextureCount() const { return tex_coords_.size(); }

  // Returns the number of vertices as an int.
  int GetVertexCount() const { return positions_.size(); }

  // Returns the number of triangles as an int.
  int GetTriangleCount() const { return triangles_.size(); }

  // Returns a const reference to the vertex positions.
  const std::vector<ion::math::Point3f>& GetPositions() const {
    return positions_;
  }

  // Returns a reference to the vertex positions.
  std::vector<ion::math::Point3f>& GetPositions() { return positions_; }

  // Returns a const reference to the i-th texture coordinate set.
  const std::vector<ion::math::Point3f>& GetTexCoords(int i) const {
    return tex_coords_.at(i);
  }

  // Returns a reference to the i-th texture coordinate set.
  std::vector<ion::math::Point3f>& GetTexCoords(int i) {
    return tex_coords_.at(i);
  }

  // Returns a const reference to the triangles.
  const std::vector<Triangle>& GetTriangles() const { return triangles_; }

  // Returns a reference to the triangles.
  std::vector<Triangle>& GetTriangles() { return triangles_; }

  // Appends a vertex to the mesh.
  void AppendVertex(const ion::math::Point3f& position,
                    absl::Span<const ion::math::Point3f> tex_coords) {
    CHECK_EQ(tex_coords_.size(), tex_coords.size());
    positions_.push_back(position);
    for (int i = 0; i < GetTextureCount(); ++i) {
      tex_coords_[i].push_back(tex_coords[i]);
    }
  }

  // Appends a triangle to the mesh.
  void AppendTriangle(const Triangle& triangle) {
    triangles_.push_back(triangle);
  }

  // Appends another |mesh| to this one. Vertices of |mesh| are appended to the
  // vertices of this mesh. Triangles are appended with the vertex indices
  // offset accordingly.
  void AppendMesh(const Mesh& mesh);

  // Transforms all vertex positions by the given matrix.
  void TransformPositions(const ion::math::Matrix4f& matrix);

 private:
  // The vertices of the mesh.
  std::vector<ion::math::Point3f> positions_;

  // Array of texture coordinates. Outer index runs over texture coordinate
  // sets, inner index runs over vertices.
  std::vector<std::vector<ion::math::Point3f>> tex_coords_;

  // The triangles of the mesh.
  std::vector<Triangle> triangles_;
};

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_MESH_H_
