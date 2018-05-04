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

// Portable implementations of FileSytem methods, based on std::fstream.

#include "seurat/base/file_system.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>

#include "absl/strings/substitute.h"

namespace seurat {
namespace base {

base::Status FileSystem::GetContents(absl::string_view path,
                                     std::string* file_contents) {
  const std::string absolute_path = GetAbsolutePath(path);

  std::ifstream local_file(absolute_path, std::ios::binary | std::ios::in);
  if (!local_file.is_open()) {
    return base::NotFoundError(
        absl::Substitute("cannot open \"$0\" for reading", absolute_path));
  }

  std::ifstream::streampos file_size = 0;
  local_file.seekg(0, local_file.end);
  file_size = local_file.tellg();
  local_file.seekg(0, local_file.beg);

  file_contents->resize(file_size);
  local_file.read(&(*file_contents)[0], file_size);

  if (!local_file.good()) {
    return base::OutOfRangeError(
        absl::Substitute("error reading from \"$0\"", absolute_path));
  }

  return base::OkStatus();
}

base::Status FileSystem::SetContents(absl::string_view path,
                                     absl::string_view file_contents) {
  const std::string absolute_path = GetAbsolutePath(path);

  std::ofstream local_file(absolute_path, std::ios::binary | std::ios::out);

  if (!local_file.is_open()) {
    return base::NotFoundError(
        absl::Substitute("cannot open \"$0\" for writing", absolute_path));
  }

  local_file.write(file_contents.data(), file_contents.size());
  if (!local_file.good()) {
    return base::OutOfRangeError(
        absl::Substitute("error writing to \"$0\"", absolute_path));
  }

  return base::OkStatus();
}

base::Status FileSystem::Match(absl::string_view pattern,
                               std::vector<std::string>* results) const {
  return base::UnimplementedError("FileSystem::Match() not implemented");
}

}  // namespace base
}  // namespace seurat
