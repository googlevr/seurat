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

#ifndef VR_SEURAT_INGEST_CLAMPING_VIEW_GROUP_LOADER_H_
#define VR_SEURAT_INGEST_CLAMPING_VIEW_GROUP_LOADER_H_

#include <memory>

#include "seurat/base/camera.h"
#include "seurat/ingest/view_group_loader.h"

namespace seurat {
namespace ingest {

// Wraps another ViewGroupLoader to clamp observed samples onto a skybox.
class ClampingViewGroupLoader : public ViewGroupLoader {
 public:
  ClampingViewGroupLoader(float skybox_radius, bool zero_is_infinite,
                          std::unique_ptr<ViewGroupLoader> delegate)
      : skybox_radius_(skybox_radius),
        zero_is_infinite_(zero_is_infinite),
        delegate_(std::move(delegate)) {}

  // ViewGroupLoader implementation.
  int GetNumViewGroups() const override {
    return delegate_->GetNumViewGroups();
  }
  base::Status LoadViewGroup(
      int view_group_index, std::vector<std::shared_ptr<base::Camera>>* cameras,
      std::vector<image::Ldi4f>* ldis) const override;

 private:
  // Wraps the given Camera to clamp ray endpoints onto the skybox.
  std::shared_ptr<base::Camera> WrapCamera(
      std::shared_ptr<base::Camera> original) const;

  // Half the side length of an origin-centered cube defining the skybox.
  const float skybox_radius_;

  // If true, depths of zero will be considered infinite.
  //
  // This is useful, for example, when processing output from a game engine
  // which cannot properly represent points at infinity in its depth buffer.
  const bool zero_is_infinite_;

  // The ViewGroupLoader to wrap.
  const std::unique_ptr<ViewGroupLoader> delegate_;
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_CLAMPING_VIEW_GROUP_LOADER_H_
