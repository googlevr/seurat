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

#include <functional>

#include "gtest/gtest.h"

namespace seurat {
namespace base {
namespace {

TEST(Status, OkStatus) {
  Status status = OkStatus();
  EXPECT_TRUE(status.ok());
}

TEST(Status, Error) {
  Status status("a");
  EXPECT_FALSE(status.ok());
  EXPECT_EQ("a", status.error_message());
}

TEST(Status, CanonicalErrors) {
  std::array<std::function<Status(absl::string_view)>, 16> factory_functions{
      {&AbortedError, &AlreadyExistsError, &CancelledError, &DataLossError,
       &DeadlineExceededError, &FailedPreconditionError, &InternalError,
       &InvalidArgumentError, &NotFoundError, &OutOfRangeError,
       &PermissionDeniedError, &UnauthenticatedError, &ResourceExhaustedError,
       &UnavailableError, &UnimplementedError, &UnknownError}};
  for (const auto& factory_function : factory_functions) {
    Status status = factory_function("test");
    EXPECT_FALSE(status.ok());
    EXPECT_EQ("test", status.error_message());
  }
}

}  // namespace
}  // namespace base
}  // namespace seurat
