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

#include "seurat/base/status_util.h"

#include "absl/strings/str_cat.h"

namespace seurat {
namespace base {

void UpdateStatus(base::Status* status, base::Status new_status) {
  if (new_status.ok()) {
    return;
  }
  if (status->ok()) {
    *status = std::move(new_status);
  } else {
    *status = base::Status(absl::StrCat(status->error_message(), "\n",
                                        new_status.error_message()));
  }
  return;
}

}  // namespace base
}  // namespace seurat
