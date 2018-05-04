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

#include "seurat/ingest/subset_view_group_loader.h"

#include "absl/strings/str_cat.h"
#include "seurat/base/camera.h"
#include "seurat/base/status.h"

namespace seurat {
namespace ingest {

int SubsetViewGroupLoader::MapViewGroupIndex(int view_group_index) const {
  float scale = static_cast<float>(delegate_->GetNumViewGroups()) /
                static_cast<float>(max_view_groups_);
  return view_group_index * scale;
}

base::Status SubsetViewGroupLoader::LoadViewGroup(
    int view_group_index, std::vector<std::shared_ptr<base::Camera>>* cameras,
    std::vector<image::Ldi4f>* ldis) const {
  if (view_group_index >= max_view_groups_) {
    return base::OutOfRangeError(
        absl::StrCat("view_group_index=", view_group_index,
                     " >= max_view_groups_=", max_view_groups_));
  }
  int mapped_index = MapViewGroupIndex(view_group_index);
  return delegate_->LoadViewGroup(mapped_index, cameras, ldis);
}

}  // namespace ingest
}  // namespace seurat
