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

#ifndef VR_SEURAT_INGEST_LDI_LOADER_TEST_UTILS_H_
#define VR_SEURAT_INGEST_LDI_LOADER_TEST_UTILS_H_

#include "seurat/image/ldi.h"
#include "seurat/ingest/ldi_loader.h"

namespace seurat {
namespace ingest {

// This is a LdiLoader class that returns a constant Ldi4f, provided at
// construction time.
class FakeLdiLoader : public LdiLoader {
 public:
  explicit FakeLdiLoader(image::Ldi4f ldi);
  ~FakeLdiLoader() override = default;

  // LdiLoader implementation.
  base::Status Load(const api::proto::Ldi& proto, base::FileSystem* file_system,
                    image::Ldi4f* ldi) const override;

 private:
  const image::Ldi4f ldi_;
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_LDI_LOADER_TEST_UTILS_H_
