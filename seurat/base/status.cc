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

#include "seurat/base/status.h"

namespace seurat {
namespace base {

Status::Status() {}

Status::Status(absl::string_view message) : error_message_(message) {}

Status OkStatus() { return Status(); }

#define STATUS_FACTORY_FUNCTION(ERROR) \
  Status ERROR(absl::string_view message) { return Status(message); }

STATUS_FACTORY_FUNCTION(AbortedError);
STATUS_FACTORY_FUNCTION(AlreadyExistsError);
STATUS_FACTORY_FUNCTION(CancelledError);
STATUS_FACTORY_FUNCTION(DataLossError);
STATUS_FACTORY_FUNCTION(DeadlineExceededError);
STATUS_FACTORY_FUNCTION(FailedPreconditionError);
STATUS_FACTORY_FUNCTION(InternalError);
STATUS_FACTORY_FUNCTION(InvalidArgumentError);
STATUS_FACTORY_FUNCTION(NotFoundError);
STATUS_FACTORY_FUNCTION(OutOfRangeError);
STATUS_FACTORY_FUNCTION(PermissionDeniedError);
STATUS_FACTORY_FUNCTION(UnauthenticatedError);
STATUS_FACTORY_FUNCTION(ResourceExhaustedError);
STATUS_FACTORY_FUNCTION(UnavailableError);
STATUS_FACTORY_FUNCTION(UnimplementedError);
STATUS_FACTORY_FUNCTION(UnknownError);

}  // namespace base
}  // namespace seurat
