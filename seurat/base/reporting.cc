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

#include <iostream>

namespace seurat {
namespace base {

// Log an info message.
void SeuratInfo(absl::string_view message) {
  std::cout << "INFO: " << message << std::endl;
}

// Log a warning message.
void SeuratWarning(absl::string_view message) {
  std::cerr << "WARNING: " << message << std::endl;
}

// Log an error message.
void SeuratError(absl::string_view message) {
  std::cerr << "ERROR: " << message << std::endl;
}

// Log a fatal error message and exit.
ABSL_ATTRIBUTE_NORETURN
void SeuratFatal(absl::string_view message) {
  std::cerr << "FATAL: " << message << std::endl;
  exit(EXIT_FAILURE);
}

}  // namespace base
}  // namespace seurat
