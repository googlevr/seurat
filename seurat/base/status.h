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

#ifndef VR_SEURAT_BASE_STATUS_H_
#define VR_SEURAT_BASE_STATUS_H_

#include <string>

#include "absl/base/attributes.h"
#include "absl/strings/string_view.h"

namespace seurat {
namespace base {

// A simple status class.
class ABSL_MUST_USE_RESULT Status {
 public:
  // Constructs a successful status with an empty (no error) message.
  Status();

  // Constructs a status with the given |message|. If the |message| string is
  // empty, the resulting status will be okay. Otherwise, it will be an error.
  explicit Status(absl::string_view message);

  ~Status() = default;

  // Compare this status to another.
  bool operator==(const Status& other) const {
    return error_message_ == other.error_message_;
  }
  bool operator!=(const Status& other) const { return !(*this == other); }

  // Returns true if this status is not an error.
  ABSL_MUST_USE_RESULT bool ok() const { return error_message_.empty(); }

  // Converts to bool.
  explicit operator bool() const { return ok(); }

  // Returns the error message.
  std::string error_message() const { return error_message_; }

 private:
  // The error message.
  std::string error_message_;
};

// Factory functions for the canonical error messages.
Status OkStatus();
Status AbortedError(absl::string_view message);
Status AlreadyExistsError(absl::string_view message);
Status CancelledError(absl::string_view message);
Status DataLossError(absl::string_view message);
Status DeadlineExceededError(absl::string_view message);
Status FailedPreconditionError(absl::string_view message);
Status InternalError(absl::string_view message);
Status InvalidArgumentError(absl::string_view message);
Status NotFoundError(absl::string_view message);
Status OutOfRangeError(absl::string_view message);
Status PermissionDeniedError(absl::string_view message);
Status UnauthenticatedError(absl::string_view message);
Status ResourceExhaustedError(absl::string_view message);
Status UnavailableError(absl::string_view message);
Status UnimplementedError(absl::string_view message);
Status UnknownError(absl::string_view message);

// Convenience macros.
#define SEURAT_RETURN_IF_ERROR(expr)        \
  if (seurat::base::Status status = expr) { \
    ;                                          \
  } else {                                     \
    return status;                             \
  }

#define SEURAT_RETURN_IF_ERROR_EXTENDED(expr, status_extension) \
  if (seurat::base::Status status = expr) {                     \
    ;                                                              \
  } else {                                                         \
    base::UpdateStatus(&status, status_extension);                 \
    return status;                                                 \
  }

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_STATUS_H_
