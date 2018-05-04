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

#ifndef VR_SEURAT_INGEST_VIEW_GROUP_LOADER_H_
#define VR_SEURAT_INGEST_VIEW_GROUP_LOADER_H_

#include <memory>
#include <vector>

#include "seurat/base/projective_camera.h"
#include "seurat/base/status.h"
#include "seurat/geometry/cube_face.h"
#include "seurat/image/ldi.h"

namespace seurat {
namespace ingest {

// Implementations of this interface provide access to structured view data,
// e.g. a JSON file describing the cameras and image files. Views are loaded in
// groups, typically the faces of a cube map or a subset of them.
class ViewGroupLoader {
 public:
  virtual ~ViewGroupLoader() = default;

  // Returns the number of view groups.
  virtual int GetNumViewGroups() const = 0;

  // Loads the |cameras| and |ldis| for the given |view_group_index|.
  //
  // If any parameters are nullptr, no value is loaded.
  virtual base::Status LoadViewGroup(
      int view_group_index, std::vector<std::shared_ptr<base::Camera>>* cameras,
      std::vector<image::Ldi4f>* ldis) const = 0;
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_VIEW_GROUP_LOADER_H_
