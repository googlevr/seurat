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

#ifndef VR_SEURAT_BASE_STATUS_UTIL_H_
#define VR_SEURAT_BASE_STATUS_UTIL_H_

#include "seurat/base/status.h"

namespace seurat {
namespace base {

// Updates a base::Status with a new Status, if the new Status is a failure
// state, in which case:
//
// * If the old Status is a success, replace it with the new Status.
// * If the old Status is a failure, then retain it but append its message with
//   the message from the new Status.
void UpdateStatus(base::Status* status, base::Status new_status);

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_STATUS_UTIL_H_
