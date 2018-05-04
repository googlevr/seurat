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

#ifndef VR_SEURAT_INGEST_VIEW_GROUP_LOADER_UTIL_H_
#define VR_SEURAT_INGEST_VIEW_GROUP_LOADER_UTIL_H_

#include "ion/math/range.h"
#include "seurat/base/status.h"
#include "seurat/ingest/view_group_loader.h"

namespace seurat {
namespace ingest {

// Computes the headbox as the bounding box around the world-space positions of
// all cameras of all view groups.
ion::math::Range3f ComputeHeadbox(const ViewGroupLoader& view_group_loader);

// Invokes |func| with each view group, in order.
//
// Additional threads are spawned to prefetch view groups in parallel with
// |func|.
base::Status ForEachViewGroupPrefetching(
    const ViewGroupLoader& loader,
    const std::function<base::Status(std::vector<std::shared_ptr<base::Camera>>,
                                     std::vector<image::Ldi4f>)>& func);

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_VIEW_GROUP_LOADER_UTIL_H_
