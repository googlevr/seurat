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

#ifndef VR_SEURAT_BASE_REPORTING_H_
#define VR_SEURAT_BASE_REPORTING_H_

#include "absl/base/port.h"
#include "absl/strings/string_view.h"

namespace seurat {
namespace base {

// Log a user visible info message.
void SeuratInfo(absl::string_view message);

// Log a user visible warning message.
void SeuratWarning(absl::string_view message);

// Log a user visible error message.
void SeuratError(absl::string_view message);

// Log a user visible fatal error message and exit.
ABSL_ATTRIBUTE_NORETURN
void SeuratFatal(absl::string_view message);

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_REPORTING_H_
