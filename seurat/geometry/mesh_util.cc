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

#include <math.h>  // For M_PI on MSVC14
#include <map>
#include <set>

#include "ion/math/vectorutils.h"
#include "absl/types/span.h"
#include "seurat/base/util.h"
#include "seurat/geometry/mesh.h"

namespace seurat {
namespace geometry {

using ion::math::Point2f;
using ion::math::Point3f;

void AppendTriangleFan(absl::Span<const Point3f> positions,
                       absl::Span<const absl::Span<const Point3f>> tex_coords,
                       Mesh* mesh) {
  CHECK_EQ(mesh->GetTextureCount(), tex_coords.size());
  if (!tex_coords.empty()) {
    for (auto tex_coords_set : tex_coords) {
      DCHECK(tex_coords_set.size() == positions.size())
          << "positions.size() = " << positions.size() << ", "
          << "tex_coords_set.size() = " << tex_coords_set.size();
    }
  }
  const int vertex_count = positions.size();
  if (vertex_count < 3) {
    return;
  }

  const int base_index = mesh->GetVertexCount();
  std::vector<Point3f> vertex_tex_coords(mesh->GetTextureCount());
  for (int i = 0; i < vertex_count; ++i) {
    for (int j = 0; j < mesh->GetTextureCount(); ++j) {
      vertex_tex_coords[j] = tex_coords[j][i];
    }
    mesh->AppendVertex(positions[i], vertex_tex_coords);
  }

  for (int i = 0; i < vertex_count - 2; ++i) {
    mesh->AppendTriangle(
        {{base_index, base_index + i + 1, base_index + i + 2}});
  }
}

float EstimateProjectedArea(const Point3f& sphere_center, const Mesh& mesh) {
  float estimated_area = 0.0f;

  absl::Span<const Point3f> positions = mesh.GetPositions();
  for (const auto& triangle : mesh.GetTriangles()) {
    std::array<Point3f, 3> vertex_positions;
    for (int i = 0; i < 3; ++i) {
      int index = triangle[i];
      vertex_positions[i] = positions[index];
    }

    for (Point3f& p : vertex_positions) {
      p = ion::math::Normalized(p - sphere_center) + Point3f::Zero();
    }

    // Area of a triangle is (1/2) * cross-product of side vectors.
    estimated_area += 0.5f * ion::math::Length(ion::math::Cross(
                                 vertex_positions[2] - vertex_positions[0],
                                 vertex_positions[1] - vertex_positions[0]));
  }

  return estimated_area / (4.0f * M_PI);
}

Mesh ToUnindexedMesh(const Mesh& indexed) {
  Mesh result(indexed.GetTextureCount());
  std::vector<Point3f> vertex_tex_coords(result.GetTextureCount());
  for (const auto& triangle : indexed.GetTriangles()) {
    int base_index = result.GetVertexCount();
    for (int i : triangle) {
      for (int j = 0; j < indexed.GetTextureCount(); ++j) {
        vertex_tex_coords[j] = indexed.GetTexCoords(j)[i];
      }
      result.AppendVertex(indexed.GetPositions()[i], vertex_tex_coords);
    }
    result.AppendTriangle({{base_index, base_index + 1, base_index + 2}});
  }
  return result;
}

std::vector<Mesh> SplitMeshByVertexCount(const Mesh& original,
                                         int max_vertex_count) {
  CHECK_GE(max_vertex_count, 3);
  int texture_count = original.GetTextureCount();
  std::vector<Mesh> meshes;
  meshes.push_back(Mesh(texture_count));
  // Maps vertices from their index in the original mesh to their assigned index
  // in the (current) new mesh.
  std::map<int, int> new_from_original_index;

  std::set<int> new_indices;

  std::vector<Point3f> vertex_tex_coords(texture_count);

  for (const Mesh::Triangle& triangle : original.GetTriangles()) {
    Mesh* cur_mesh = &meshes.back();

    // Unique new indices.
    new_indices.clear();
    for (int index : triangle) {
      new_indices.insert(index);
    }

    int vertices_to_add = 0;
    for (int index : new_indices) {
      if (new_from_original_index.count(index) == 0) {
        vertices_to_add++;
      }
    }

    // If the vertex count would overflow, create another mesh.
    if (vertices_to_add + cur_mesh->GetVertexCount() > max_vertex_count) {
      new_from_original_index.clear();
      meshes.push_back(Mesh(texture_count));
      cur_mesh = &meshes.back();
    }

    for (int index : triangle) {
      if (new_from_original_index.count(index) == 0) {
        int new_index = cur_mesh->GetVertexCount();
        for (int j = 0; j < texture_count; ++j) {
          vertex_tex_coords[j] = original.GetTexCoords(j)[index];
        }
        cur_mesh->AppendVertex(original.GetPositions()[index],
                               vertex_tex_coords);
        new_from_original_index[index] = new_index;
      }
    }

    cur_mesh->AppendTriangle({{new_from_original_index[triangle[0]],
                               new_from_original_index[triangle[1]],
                               new_from_original_index[triangle[2]]}});
  }
  return meshes;
}

}  // namespace geometry
}  // namespace seurat
