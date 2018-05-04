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

#include "seurat/ingest/ldi_loader_test_utils.h"

#include "seurat/base/status.h"

namespace seurat {
namespace ingest {

FakeLdiLoader::FakeLdiLoader(image::Ldi4f ldi) : ldi_(std::move(ldi)) {}

base::Status FakeLdiLoader::Load(const api::proto::Ldi& proto,
                                 base::FileSystem* file_system,
                                 image::Ldi4f* ldi) const {
  *ldi = ldi_;
  return base::OkStatus();
}

}  // namespace ingest
}  // namespace seurat
