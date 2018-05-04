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

#ifndef VR_SEURAT_INGEST_PROTO_VIEW_GROUP_LOADER_H_
#define VR_SEURAT_INGEST_PROTO_VIEW_GROUP_LOADER_H_

#include <memory>
#include <vector>

#include "seurat/api/api.pb.h"
#include "seurat/base/camera.h"
#include "seurat/base/file_system.h"
#include "seurat/geometry/cube_face.h"
#include "seurat/image/ldi.h"
#include "seurat/ingest/ldi_loader.h"
#include "seurat/ingest/view_group_loader.h"

namespace seurat {
namespace ingest {

class ProtoViewGroupLoader : public ViewGroupLoader {
 public:
  ProtoViewGroupLoader(const api::proto::Capture& capture, int thread_count,
                       std::unique_ptr<LdiLoader> ldi_loader,
                       std::shared_ptr<base::FileSystem> file_system);
  ~ProtoViewGroupLoader() override = default;

  // ViewGroupLoader implementation.
  int GetNumViewGroups() const override;
  base::Status LoadViewGroup(
      int view_group_index, std::vector<std::shared_ptr<base::Camera>>* cameras,
      std::vector<image::Ldi4f>* ldis) const override;

 private:
  // The proto defining the view groups.
  const api::proto::Capture capture_;

  // Controls the number of threads the loader runs to load views.
  const int thread_count_;

  // The LDI loader.
  std::unique_ptr<LdiLoader> ldi_loader_;

  // The file system used to load LDIs.
  std::shared_ptr<base::FileSystem> file_system_;
};

}  // namespace ingest
}  // namespace seurat

#endif  //  VR_SEURAT_INGEST_PROTO_VIEW_GROUP_LOADER_H_
