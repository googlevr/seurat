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

#ifndef VR_SEURAT_GEOMETRY_MESH_OBJ_IO_H_
#define VR_SEURAT_GEOMETRY_MESH_OBJ_IO_H_

#include <string>

#include "seurat/base/file_system.h"
#include "seurat/base/status.h"
#include "seurat/geometry/mesh.h"

namespace seurat {
namespace geometry {

class MeshObjIo {
 public:
  // Writes a |mesh| as OBJ to the string |out_filedata|.
  static void WriteObj(std::string* out_filedata, const Mesh& mesh);

  // Writes a |mesh| as OBJ to the file with name |filename|. Returns the
  // status of the write operation.
  static base::Status WriteObjToFile(base::FileSystem* file_system,
                                     const std::string& filename,
                                     const Mesh& mesh);
};

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_MESH_OBJ_IO_H_
