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

#include "seurat/ingest/tone_mapping_view_group_loader.h"

#include "ion/math/vector.h"
#include "seurat/base/camera.h"
#include "seurat/base/status.h"
#include "seurat/image/ldi.h"

namespace seurat {
namespace ingest {

using base::Camera;
using image::Ldi4f;
using ion::math::Point3f;
using ion::math::Vector3f;

base::Status ToneMappingViewGroupLoader::LoadViewGroup(
    int view_group_index, std::vector<std::shared_ptr<Camera>>* cameras,
    std::vector<Ldi4f>* ldis) const {
  base::Status status =
      delegate_->LoadViewGroup(view_group_index, cameras, ldis);

  if (ldis != nullptr) {
    for (Ldi4f& ldi : *ldis) {
      tone_mapper_->ProcessColors(ldi.GetMutableColors());
    }
  }

  return status;
}

}  // namespace ingest
}  // namespace seurat
