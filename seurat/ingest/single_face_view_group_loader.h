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

#ifndef VR_SEURAT_INGEST_SINGLE_FACE_VIEW_GROUP_LOADER_H_
#define VR_SEURAT_INGEST_SINGLE_FACE_VIEW_GROUP_LOADER_H_

#include <memory>

#include "seurat/base/camera.h"
#include "seurat/ingest/view_group_loader.h"

namespace seurat {
namespace ingest {

// Wraps another ViewGroupLoader to restrict samples to only those which project
// onto the specified cube face.  All other samples are culled.
class SingleFaceViewGroupLoader : public ViewGroupLoader {
 public:
  SingleFaceViewGroupLoader(const geometry::CubeFace cube_face,
                            int thread_count,
                            std::unique_ptr<ViewGroupLoader> delegate)
      : cube_face_(cube_face),
        thread_count_(thread_count),
        delegate_(std::move(delegate)) {}

  // ViewGroupLoader implementation.
  int GetNumViewGroups() const override {
    return delegate_->GetNumViewGroups();
  }
  base::Status LoadViewGroup(
      int view_group_index, std::vector<std::shared_ptr<base::Camera>>* cameras,
      std::vector<image::Ldi4f>* ldis) const override;

 private:
  // The cube face to clip points against.
  const geometry::CubeFace cube_face_;

  // The number of threads to use.
  const int thread_count_;

  // The ViewGroupLoader to wrap.
  const std::unique_ptr<ViewGroupLoader> delegate_;
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_SINGLE_FACE_VIEW_GROUP_LOADER_H_
