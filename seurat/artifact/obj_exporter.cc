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

#include "seurat/artifact/obj_exporter.h"

#include "absl/strings/str_cat.h"
#include "seurat/base/status.h"
#include "seurat/geometry/mesh.h"
#include "seurat/geometry/mesh_obj_io.h"

namespace seurat {
namespace artifact {

using geometry::Mesh;
using geometry::MeshObjIo;

base::Status ObjExporter::Process(Artifact* artifact) const {
  if (!artifact->mesh) {
    return base::UnimplementedError("No mesh to export.");
  }
  std::string filename = absl::StrCat(basename_, ".obj");
  return MeshObjIo::WriteObjToFile(filesystem_.get(), filename,
                                   *artifact->mesh);
}

}  // namespace artifact
}  // namespace seurat
