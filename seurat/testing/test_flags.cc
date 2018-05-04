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

#include "seurat/testing/test_flags.h"

#include <stdlib.h>
#include <string>

#include "ion/base/logging.h"

namespace seurat {
namespace testing {

std::string GetTestSrcdir() {
  const char* value = getenv("TEST_SRCDIR");
  if (value == nullptr) {
    LOG(FATAL) << "Tests referencing src files should be run via 'bazel test'";
  }
  return std::string(value);
}

std::string GetTestTmpdir() {
  const char* value = getenv("TEST_TMPDIR");
  if (value == nullptr) {
#if defined(OS_WINDOWS)
    LOG(FATAL)
        << "Tests requiring a tmp directory should be run via 'bazel test'";
#else
    // On linux, just default to /tmp for convenience.
    return "/tmp";
#endif
  }
  return std::string(value);
}

}  // namespace testing
}  // namespace seurat
