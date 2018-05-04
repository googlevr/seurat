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

#include "seurat/geometry/mesh_obj_io.h"

#include <string>

#include "seurat/base/file_system.h"
#include "seurat/base/reporting.h"
#include "seurat/geometry/mesh.h"

namespace seurat {
namespace geometry {

using ion::math::Point2f;
using ion::math::Point3f;

namespace {

void WriteObjInner(std::string* out_filedata, const Mesh& mesh) {
  // Write vertex positions
  for (const Point3f& position : mesh.GetPositions()) {
    *out_filedata += "v " +                               //
                     std::to_string(position[0]) + " " +  //
                     std::to_string(position[1]) + " " +  //
                     std::to_string(position[2]) + "\n";  //
  }

  // Write texture coordinates
  if (mesh.GetTextureCount() > 0) {
    bool has_homogeneous_texture_coordinates = false;
    for (const Point3f& tex_coord : mesh.GetTexCoords(0)) {
      if (tex_coord[2] != 1.0f) {
        has_homogeneous_texture_coordinates = true;
        *out_filedata += "vt " +                               //
                         std::to_string(tex_coord[0]) + " " +  //
                         std::to_string(tex_coord[1]) + " " +  //
                         std::to_string(tex_coord[2]) + "\n";  //
      } else {
        *out_filedata += "vt " +                               //
                         std::to_string(tex_coord[0]) + " " +  //
                         std::to_string(tex_coord[1]) + "\n";  //
      }
    }
    if (has_homogeneous_texture_coordinates) {
      base::SeuratWarning(
          "Exporting OBJ with homogeneous (projective) texture coordinates. "
          "Many applications don't support this format properly.");
    }
  }

  // Write triangles
  for (const Mesh::Triangle& triangle : mesh.GetTriangles()) {
    // OBJ uses 1-based indices.
    const int i0 = triangle[0] + 1;
    const int i1 = triangle[1] + 1;
    const int i2 = triangle[2] + 1;
    *out_filedata += "f " +                                                 //
                     std::to_string(i0) + "/" + std::to_string(i0) + " " +  //
                     std::to_string(i1) + "/" + std::to_string(i1) + " " +  //
                     std::to_string(i2) + "/" + std::to_string(i2) + "\n";  //
  }
}

}  // namespace

void MeshObjIo::WriteObj(std::string* out_filedata, const Mesh& mesh) {
  if (mesh.GetTextureCount() > 1) {
    LOG(ERROR) << "Mesh exported to OBJ file format has more than one texture "
                  "coordinate set. Writing out only the first set.";
  }
  WriteObjInner(out_filedata, mesh);
}

base::Status MeshObjIo::WriteObjToFile(base::FileSystem* file_system,
                                       const std::string& filename,
                                       const Mesh& mesh) {
  CHECK(!filename.empty()) << "filename must not be empty";
  std::string out_obj_data;
  if (mesh.GetTextureCount() > 1) {
    LOG(ERROR) << "Mesh exported to .obj file " << filename
               << " has more than one texture coordinate set. Writing out only "
                  "the first set.";
  }
  WriteObjInner(&out_obj_data, mesh);
  return file_system->SetContents(filename, out_obj_data);
}

}  // namespace geometry
}  // namespace seurat
