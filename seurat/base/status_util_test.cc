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

#include "gtest/gtest.h"
#include "absl/strings/str_cat.h"

namespace seurat {
namespace base {
namespace {

TEST(StatusUtilTest, TestUpdateStatus) {
  base::Status status_cases[][3] = {
      // OK + OK == OK
      {base::OkStatus(), base::OkStatus(), base::OkStatus()},
      // OK + err2 == err2
      {base::OkStatus(), base::Status("foo"), base::Status("foo")},
      // err1 + OK == err1
      {base::Status("foo"), base::OkStatus(), base::Status("foo")},
      // err1 + err2 == err1 + err2.msg
      {base::Status("foo"), base::Status("bar"), base::Status("foo\nbar")},
  };
  for (const auto& row : status_cases) {
    base::Status actual = row[0];
    const base::Status new_status = row[1];
    const base::Status expected = row[2];
    UpdateStatus(&actual, new_status);
    EXPECT_EQ(expected, actual);
  }
}

}  // namespace
}  // namespace base
}  // namespace seurat
