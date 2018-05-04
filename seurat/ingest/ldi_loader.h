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

#ifndef VR_SEURAT_INGEST_LDI_LOADER_H_
#define VR_SEURAT_INGEST_LDI_LOADER_H_

#include "seurat/api/image.pb.h"
#include "seurat/base/file_system.h"
#include "seurat/base/status.h"
#include "seurat/image/ldi.h"

namespace seurat {
namespace ingest {

class LdiLoader {
 public:
  LdiLoader() = default;
  virtual ~LdiLoader() = default;

  // Load an Ldi4f from a file on a FileSystem, using a api::proto::Ldi
  // specification.
  virtual base::Status Load(const api::proto::Ldi& proto,
                            base::FileSystem* file_system,
                            image::Ldi4f* ldi) const;
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_LDI_LOADER_H_
