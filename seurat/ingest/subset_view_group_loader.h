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

#ifndef VR_SEURAT_INGEST_SUBSET_VIEW_GROUP_LOADER_H_
#define VR_SEURAT_INGEST_SUBSET_VIEW_GROUP_LOADER_H_

#include <memory>

#include "seurat/ingest/view_group_loader.h"

namespace seurat {
namespace ingest {

// Wraps another ViewGroupLoader to select a subset of all views, up to a
// specified maximum number.
class SubsetViewGroupLoader : public ViewGroupLoader {
 public:
  SubsetViewGroupLoader(std::unique_ptr<ViewGroupLoader> delegate,
                        int max_view_groups)
      : delegate_(std::move(delegate)),
        max_view_groups_(
            std::min(max_view_groups, delegate_->GetNumViewGroups())) {
    CHECK_GE(max_view_groups_, 0);
  }

  // ViewGroupLoader implementation.
  int GetNumViewGroups() const override { return max_view_groups_; }
  base::Status LoadViewGroup(
      int view_group_index, std::vector<std::shared_ptr<base::Camera>>* cameras,
      std::vector<image::Ldi4f>* ldis) const override;

 private:
  // Maps the |view_group_index| to the range of the delegate.
  int MapViewGroupIndex(int view_gropup_index) const;

  // The ViewGroupLoader to wrap.
  const std::unique_ptr<ViewGroupLoader> delegate_;

  // The number of view groups returned by this loader.
  const int max_view_groups_;
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_SUBSET_VIEW_GROUP_LOADER_H_
