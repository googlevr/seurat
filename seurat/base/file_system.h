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

#ifndef VR_SEURAT_BASE_FILE_SYSTEM_H_
#define VR_SEURAT_BASE_FILE_SYSTEM_H_

#include <string>
#include <utility>
#include <vector>

#include "absl/strings/string_view.h"
#include "seurat/base/status.h"

namespace seurat {
namespace base {

// Encapsulate access to various host file systems for Seurat.
//
// FileSystem supports both Windows and Unix platforms, which have different
// path conventions.  Cross-compatibility support is as follows:
//
// * FileSystem on Windows supports native Windows paths, as well as Unix paths.
//   JoinPath() will successfully join Windows and Unix paths; when the result
//   is cleaned with CleanPath(), a "clean" Windows path is returned.
//   * Unix relative paths are converted to Windows relative paths.
//   * Unix absolute paths are resolved to root at the drive letter of current
//     working directory, in keeping with Windows convention.  For example, with
//     a current working directory of "D:\foo", "/bar" resolves to "D:\bar".
// * FileSystem on Unix supports Unix paths.
//
// Thus, paths that are expected to interchange between platforms should be
// written in Unix format.
class FileSystem {
 public:
  // Construct a file system capable of interacting with the machine's local
  // file system.  |root_path| is the root for relative paths used by this
  // FileSystem; all such paths are interpreted relative to |root_path|.
  // Absolute paths are still resolved absolutely.
  explicit FileSystem(absl::string_view root_path);

  // Disallow copy and assign, but allow move.
  FileSystem(const FileSystem& other) = delete;
  FileSystem(FileSystem&& other) = default;
  FileSystem& operator=(const FileSystem& other) = delete;
  FileSystem& operator=(FileSystem&& other) = default;
  ~FileSystem() = default;

  // Returns |path| as an absolute path.  Relative paths are resolved relative
  // to this FileSystem's root path; absolute paths are resolved absolutely.
  std::string GetAbsolutePath(absl::string_view path) const;

  // Attempt to read the entire file from source location |path| into
  // |file_contents| and return operation status.
  base::Status GetContents(absl::string_view path, std::string* file_contents);

  // Attempt to write |file_contents| to a file at destination location
  // |path| and return operation status.
  base::Status SetContents(absl::string_view path,
                           absl::string_view file_contents);

  // Expands globs in the given pattern, appending results to the given vector.
  //
  // Matched results are absolute paths.
  base::Status Match(absl::string_view pattern,
                     std::vector<std::string>* results) const;

  // -- Portable file:: Path Operations --

  // Returns |path|, "cleaned" by collapsing redundant path separators,
  // resolving any "." or ".." components, and removing any trailing path
  // separators.  |path| may be a relative or absolute path.
  // See: file::CleanPath().
  static std::string CleanPath(absl::string_view path);

  // Separate |path| into two string fragments, the dirname and the
  // basename, around the last path directory separator.
  // See: file::SplitPath().
  static std::pair<absl::string_view, absl::string_view> SplitPath(
      absl::string_view path);

  // Separate |path|'s file base name around the extension delimiter.
  // See: file::SplitBasename().
  static std::pair<absl::string_view, absl::string_view> SplitBasename(
      absl::string_view path);

  // Join path fragments together with the appropriate directory separator.
  // Relative paths are joined together, but absolute paths cause all previous
  // fragments to be ignored.
  // See: file::JoinPathRespectAbsolute().
  template <typename... Args>
  static std::string JoinPath(Args&&... args);

 private:
  // Base case for the recursive variadic template expansion.
  static void JoinPathInternal(std::vector<std::string::value_type>* path,
                               absl::string_view part);

  // Variadic template version for joining an arbitrary number of path
  // components.
  template <typename... Args>
  static void JoinPathInternal(std::vector<std::string::value_type>* path,
                               absl::string_view first, Args&&... args) {
    JoinPathInternal(path, first);
    JoinPathInternal(path, std::forward<Args>(args)...);
  }

  // The normal root path of this FileSystem instance.
  std::string root_path_;
};

template <typename... Args>
// static
std::string FileSystem::JoinPath(Args&&... args) {
  std::vector<std::string::value_type> path;
  JoinPathInternal(&path, std::forward<Args>(args)...);
  return std::string(path.begin(), path.end());
}

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_FILE_SYSTEM_H_
