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

#include "seurat/base/file_system.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>

#include "ion/base/logging.h"
#include "absl/strings/substitute.h"

#if defined(ION_PLATFORM_WINDOWS)
// Necessary to do some Windows path manipulations.
#include <Windows.h>
#include <cctype>
#endif

namespace seurat {
namespace base {
namespace {

constexpr std::string::value_type kPathSeparator =
#if defined(ION_PLATFORM_WINDOWS)
    '\\';
#else
    '/';
#endif

// Returns true iff |c| is a valid path separator.
bool IsPathSeparator(std::string::value_type c) {
  return
#if defined(ION_PLATFORM_WINDOWS)
      // Also support the Unix-style path separator on Windows.
      (c == '/') ||
#endif
      (c == kPathSeparator);
}

// Returns true iff |path| (all of it) represents a valid filesystem root.
bool IsRootPath(absl::string_view path) {
  if (path.empty()) {
    return false;
  }

  if (std::all_of(path.begin(), path.end(), IsPathSeparator)) {
    // Unix-style absolute path.
    return true;
  }

#if defined(ION_PLATFORM_WINDOWS)
  auto iter = path.begin();
  // Also support Windows-style "C:\" absolute paths.
  if (iter == path.end() || !std::isalpha(*iter)) {
    return false;
  }
  ++iter;
  if (iter == path.end() || *iter != ':') {
    return false;
  }
  ++iter;
  if (iter != path.end() && std::all_of(iter, path.end(), IsPathSeparator)) {
    return true;
  }
#endif

  return false;
}

// Returns true iff |path| is an absolute path.
bool IsAbsolutePath(absl::string_view path) {
  if (path.empty()) {
    return false;
  }
  if (IsPathSeparator(path.front())) {
    // Unix-style absolute path.
    return true;
  }

#if defined(ION_PLATFORM_WINDOWS)
  if (path.size() >= 3 && std::isalpha(path[0]) && path[1] == ':' &&
      IsPathSeparator(path[2])) {
    // Windows-style absolute path.
    return true;
  }
#endif

  return false;
}

}  // namespace

FileSystem::FileSystem(absl::string_view root_path) : root_path_(root_path) {}

std::string FileSystem::GetAbsolutePath(absl::string_view path) const {
  return CleanPath(JoinPath(root_path_, path));
}

// static
std::string FileSystem::CleanPath(absl::string_view path) {
  // Clean |path| by:
  //
  // * Collapsing redundant path separators.
  // * Resolving . and .. components.
  // * Removing any trailing path separators.
  std::string clean_path;
  clean_path.reserve(path.size());

  auto iter = path.begin();
#if defined(ION_PLATFORM_WINDOWS)
  if (iter != path.end() && IsPathSeparator(*iter)) {
    // On Windows, a root path like '/' is resolved to the drive letter of the
    // current working directory.  So we prepend this drive letter to
    // |clean_path|.
    clean_path.reserve(path.size() + 2);
    WCHAR pwd[MAX_PATH];
    ::GetCurrentDirectoryW(MAX_PATH, pwd);
    const std::string::value_type drive_letter = std::toupper(pwd[0]);
    clean_path.push_back(drive_letter);
    clean_path.push_back(':');
    clean_path.push_back('\\');
    ++iter;
    LOG(WARNING) << "FileSystem::CleanPath(): using drive root \""
                 << drive_letter << ":\\\" for path \"" << path << "\"";
  }
#endif

  while (iter != path.end()) {
    if (*iter == '.') {
      auto dot_iter = iter;
      ++dot_iter;
      if (dot_iter == path.end() || IsPathSeparator(*dot_iter)) {
        // This is a "." component.  Skip it and any trailing path separators.
        while (dot_iter != path.end() && IsPathSeparator(*dot_iter)) {
          ++dot_iter;
        }
        iter = dot_iter;
        continue;
      } else if (*dot_iter == '.') {
        ++dot_iter;
        if (dot_iter == path.end() || IsPathSeparator(*dot_iter)) {
          // This is a ".." component.  Skip it and any trailing path
          // separators.
          while (dot_iter != path.end() && IsPathSeparator(*dot_iter)) {
            ++dot_iter;
          }
          iter = dot_iter;

          // We need to go up one level in |clean_path|.  If |clean_path| is the
          // root, it stays the root.  If it is empty or it ends in a "..", then
          // we just append another "/.." to it; otherwise we drop the last path
          // component.
          if (!clean_path.empty()) {
            if (IsRootPath(clean_path)) {
              // This is the root, it stays the root.
              continue;
            }
            // Find the last path component in |clean_path|.
            auto norm_riter = clean_path.rbegin();
            ++norm_riter;
            while (norm_riter != clean_path.rend() &&
                   !IsPathSeparator(*norm_riter)) {
              ++norm_riter;
            }

            const size_t norm_pos =
                std::distance(clean_path.begin(), norm_riter.base());
            if (clean_path.compare(norm_pos, clean_path.size() - norm_pos,
                                   "..") != 0) {
              // The last path component was not "..", so we erase it.
              clean_path.erase(norm_riter.base(), clean_path.end());
              continue;
            } else {
              // The last path component was "..", so append another path
              // separator  before appending ".." below.
              clean_path.push_back(kPathSeparator);
            }
          }
          // |clean_path| was empty, or ended in a "..", so we just append a
          // "../".
          clean_path.push_back('.');
          clean_path.push_back('.');
          clean_path.push_back(kPathSeparator);
          continue;
        }
      }
    }
    // Either the component did not start in a '.', or it did but is not a "."
    // or "..", so it is {0 or more non-separator} + {0 or more separator}.
    // Collapse this to {non-separator sequence} + {0 or one separator}.
    while (!IsPathSeparator(*iter)) {
      // Append the {non-separator sequence}
      clean_path.push_back(*iter);
      ++iter;
      if (iter == path.end()) {
        break;
      }
    }
    if (iter != path.end()) {
      // There is at least one separator.
      clean_path.push_back(kPathSeparator);
      ++iter;
      while (iter != path.end() && IsPathSeparator(*iter)) {
        ++iter;
      }
    }
  }

  if (clean_path.empty()) {
    // Replace the empty path with '.'.
    clean_path.push_back('.');
  } else if (!IsRootPath(clean_path) && IsPathSeparator(clean_path.back())) {
    // If this is not a root path, remove a trailing separator.
    clean_path.pop_back();
  }
  return clean_path;
}

// static
std::pair<absl::string_view, absl::string_view> FileSystem::SplitPath(
    absl::string_view path) {
  auto iter = path.rbegin();
  while (iter != path.rend() && !IsPathSeparator(*iter)) {
    ++iter;
  }
  if (iter == path.rend()) {
    // No separator found.
    return {{}, path};
  }
  if (IsRootPath(absl::string_view(path.begin(),
                                   std::distance(path.begin(), iter.base())))) {
    // Separator found, but it's part of the root path.  The separator should be
    // part of the dirname.
    return {
        absl::string_view(path.begin(),
                          std::distance(path.begin(), iter.base())),
        absl::string_view(iter.base(), std::distance(iter.base(), path.end()))};
  }
  // Separator found, and it's not part of the root path.  The separator should
  // not be part of the dirname.
  return {
      absl::string_view(path.begin(),
                        std::distance(path.begin(), std::prev(iter.base()))),
      absl::string_view(iter.base(), std::distance(iter.base(), path.end()))};
}

// static
std::pair<absl::string_view, absl::string_view> FileSystem::SplitBasename(
    absl::string_view path) {
  path = SplitPath(path).second;

  const auto pos = path.find_last_of('.');
  if (pos == absl::string_view::npos)
    return std::make_pair(path, absl::ClippedSubstr(path, path.size(), 0));
  return std::make_pair(path.substr(0, pos),
                        absl::ClippedSubstr(path, pos + 1));
}

// static
void FileSystem::JoinPathInternal(std::vector<std::string::value_type>* path,
                                  absl::string_view part) {
  if (IsAbsolutePath(part)) {
    path->assign(part.begin(), part.end());
  } else {
    if (!path->empty() && !IsPathSeparator(path->back())) {
      path->push_back(kPathSeparator);
    }
    path->insert(path->end(), part.begin(), part.end());
  }
}

}  // namespace base
}  // namespace seurat
