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

#include "ion/math/matrix.h"
#include "ion/math/transformutils.h"

namespace seurat {
namespace geometry {

using ion::math::Matrix3f;
using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point3f;

void Mesh::AppendMesh(const Mesh& other) {
  CHECK_EQ(GetTextureCount(), other.GetTextureCount());
  int offset = GetVertexCount();
  positions_.insert(positions_.end(), other.positions_.begin(),
                    other.positions_.end());
  for (int i = 0; i < GetTextureCount(); ++i) {
    tex_coords_[i].insert(tex_coords_[i].end(), other.tex_coords_[i].begin(),
                          other.tex_coords_[i].end());
  }
  for (const Triangle& triangle : other.GetTriangles()) {
    this->AppendTriangle(Triangle{
        {triangle[0] + offset, triangle[1] + offset, triangle[2] + offset}});
  }
}

void Mesh::TransformPositions(const Matrix4f& matrix) {
  for (Point3f& position : positions_) {
    position = ion::math::ProjectPoint(matrix, position);
  }
}

}  // namespace geometry
}  // namespace seurat
