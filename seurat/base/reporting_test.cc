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

#include "seurat/base/reporting.h"

#include "gtest/gtest.h"

namespace seurat {
namespace base {
namespace {

TEST(Logging, Info) { SeuratInfo("Test Info"); }
TEST(Logging, Warning) { SeuratWarning("Test Warning"); }
TEST(Logging, Error) { SeuratError("Test Error"); }
TEST(Logging, Fatal) {
  EXPECT_EXIT(SeuratFatal("Test Fatal"),
              ::testing::ExitedWithCode(EXIT_FAILURE), "FATAL: Test Fatal");
}

}  // namespace
}  // namespace base
}  // namespace seurat
