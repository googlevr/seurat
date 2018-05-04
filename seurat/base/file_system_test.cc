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
#include <string>

#include "ion/port/fileutils.h"
#include "gtest/gtest.h"
#include "seurat/base/status.h"
#include "seurat/testing/test_flags.h"

using seurat::base::FileSystem;

namespace seurat {
namespace base {
namespace {

#if defined(ION_PLATFORM_WINDOWS)
std::string PrependCurrentDriveLetter(std::string s) {
  const std::string cwd = ion::port::GetCurrentWorkingDirectory();
  s.insert(s.begin(), cwd[0]);
  return s;
}
#endif

// Test GetAbsolutePath.
TEST(FileSystem, GetAbsolutePath) {
  struct GetAbsolutePathTestCase {
    std::string input;
    std::string output;
  };

  const GetAbsolutePathTestCase kGetAbsolutePathTests[] = {
#if defined(ION_PLATFORM_WINDOWS)
    {{""}, {"C:\\foo\\bar"}},           // Empty
    {{"."}, {"C:\\foo\\bar"}},          // Current directory
    {{"./"}, {"C:\\foo\\bar"}},         // Current directory, trailing '\\'
    {{".."}, {"C:\\foo"}},              // Up-level
    {{"../"}, {"C:\\foo"}},             // Up-level, trailing '\\'
    {{"../../"}, {"C:\\"}},             // Up-level to the drive
    {{"baz"}, {"C:\\foo\\bar\\baz"}},   // Relative path
    {{"baz/"}, {"C:\\foo\\bar\\baz"}},  // Relative path, trailing '\\'
    {{"/"}, PrependCurrentDriveLetter(":\\")},  // Gets the drive letter.
    {{"/d/baz"}, PrependCurrentDriveLetter(":\\d\\baz")},  // Absolute path
#else
    {{""}, {"/foo/bar"}},          // Empty
    {{"."}, {"/foo/bar"}},         // Current directory
    {{"./"}, {"/foo/bar"}},        // Current directory, trailing '/'
    {{".."}, {"/foo"}},            // Up-level
    {{"../"}, {"/foo"}},           // Up-level, trailing '/'
    {{"../../"}, {"/"}},           // Up-level past the root
    {{"baz"}, {"/foo/bar/baz"}},   // Relative path
    {{"baz/"}, {"/foo/bar/baz"}},  // Relative path, trailing '/'
    {{"/"}, {"/"}},                // Root
    {{"/baz/qux"}, {"/baz/qux"}},  // Absolute path
#endif
  };

#if defined(ION_PLATFORM_WINDOWS)
  FileSystem file_system("C:\\foo\\bar");
#else
  FileSystem file_system("/foo/bar");
#endif

  for (auto& test_case : kGetAbsolutePathTests) {
    auto absolute_result = file_system.GetAbsolutePath(test_case.input);
    EXPECT_EQ(test_case.output, absolute_result) << test_case.input;
  }
}

// Test file save failure.
TEST(FileSystem, SetContentsFailure) {
  FileSystem files(seurat::testing::GetTestTmpdir());
  base::Status result;
  std::string file_data("This data won't come back.");
  // Note: since this will be happening in the temp folder, given to us by the
  // test system, we should expect that it is writable, but the sub folder path
  // indicated here won't actually exist, so opening the file will fail.
  result = files.SetContents(
      "dir_does_not_exist/temp_file_that_cannot_be_written.txt", file_data);
  EXPECT_FALSE(result.ok());

  // The file system shouldn't give any data back for our failing location!
  std::string save_readback;
  result = files.GetContents(
      "dir_does_not_exist/temp_file_that_cannot_be_written.txt",
      &save_readback);
  EXPECT_FALSE(result.ok());
  EXPECT_TRUE(save_readback.empty());
}

// Test file open failure.
TEST(FileSystem, GetContentsFailure) {
  FileSystem files("/");
  base::Status result;
  std::string dest;
  result = files.GetContents("temp_file_that_does_not_exist.txt", &dest);
  EXPECT_FALSE(result.ok());
  EXPECT_TRUE(dest.empty());
}

// Test file save/read.
TEST(FileSystem, SetGetContents) {
  FileSystem files(seurat::testing::GetTestTmpdir());
  std::string tmp_file("test_save_read.txt");

  const std::string source_data{
      "File contents that must pass through SetContents to GetContents"};
  base::Status save_result = files.SetContents(tmp_file, source_data);
  EXPECT_TRUE(save_result.ok());

  std::string dest_data;
  base::Status read_result = files.GetContents(tmp_file, &dest_data);
  EXPECT_TRUE(read_result.ok());

  EXPECT_EQ(source_data, dest_data);
}

TEST(FileSystem, CleanPath) {
  struct CleanTestCase {
    std::string input;
    std::string output;
  };

  const CleanTestCase kCleanTests[] = {
    {{""}, {"."}},                               // Empty -> "."
    {{"./foo/bar"}, {"foo/bar"}},                // Relative path
    {{"foo/.bar/baz"}, {"foo/.bar/baz"}},        // Dot file
    {{"foo/bar.baz/qux"}, {"foo/bar.baz/qux"}},  // Dot file
    {{"foo/..bar/baz"}, {"foo/..bar/baz"}},      // Double dot file
#if defined(ION_PLATFORM_WINDOWS)
    {{"foo/bar"}, {"foo\\bar"}},  // Windows: convert Unix-style
    {{"/foo/bar"}, PrependCurrentDriveLetter(":\\foo\\bar")},  // Absolute
    {{"/foo/../.."}, PrependCurrentDriveLetter(":\\")},  // Up-level past root
    {{"C:\\foo\\.."}, {"C:\\"}},                         // Windows: to the root
    {{"C:\\foo\\..\\.."}, {"C:\\"}},  // Windows: up-level past the root
    {{"foo\\bar\\.bashrc"},
     {"foo\\bar\\.bashrc"}},  // Dots as part of filename. N.B. using Bash
                              // style volume labels as one might find working
                              // in a Git repo on Windows.
    {{"Z:\\foo\\bar\\..unusual"},
     {"Z:\\foo\\bar\\..unusual"}},  // Dots as part of filename.
    {{"Z:\\foo\\bar\\unusual.."},
     {"Z:\\foo\\bar\\unusual.."}},  // Dots as part of filename.
    {{"C:\\foo\\Assets\\..\\sub\\data.txt"},
     {"C:\\foo\\sub\\data.txt"}},  // Windows: Interior .. with trailing path.
    {{"C:\\foo/bar\\\\baz//qux/"}, {"C:\\foo\\bar\\baz\\qux"}},  // Everything
#else
    {{"/"}, {"/"}},                       // Root directory
    {{"/foo/."}, {"/foo"}},               // Current directory
    {{"/foo/./"}, {"/foo"}},              // Current directory, trailing '/'
    {{"////foo///bar//"}, {"/foo/bar"}},  // Collapsing '/'
    {{"/foo/.."}, {"/"}},                 // Relative traverse to root.
    {{"/foo/../"}, {"/"}},                // Relative traverse to root.
    {{"/foo/../file.txt"},
     {"/file.txt"}},  // Relative traverse to root with file.
    {{"/foo/bar/.bashrc"}, {"/foo/bar/.bashrc"}},  // Dots as part of filename.
    {{"/foo/bar/..unusual"},
     {"/foo/bar/..unusual"}},  // Dots as part of filename.
    {{"/foo/bar/unusual.."},
     {"/foo/bar/unusual.."}},  // Dots as part of filename.
    {{"/foo/bar/un..usual"},
     {"/foo/bar/un..usual"}},            // Dots as part of filename.
    {{"/foo/bar/.."}, {"/foo"}},         // Up-level
    {{"/foo/bar/../"}, {"/foo"}},        // Up-level, trailing '/'
    {{"/foo/bar/../../.."}, {"/"}},      // Up-level past the root.
    {{"/foo/bar/../../../../"}, {"/"}},  // Up-level further past the root.
    {{"/foo/bar/../../bar"}, {"/bar"}},  // Up-level past the root and back.
    {{"/foo/bar/../sub/final.txt"},
     {"/foo/sub/final.txt"}},         // Interior .. with trailing path.
    {{"foo/../../bar"}, {"../bar"}},  // Up-level past relative path.
    {{"foo/bar/..///./.././../baz/./qux"}, {"../baz/qux"}},  // Everything!
#endif
  };

  for (auto& test_case : kCleanTests) {
    auto clean_result = FileSystem::CleanPath(test_case.input);
    std::string test_case_output = test_case.output;
#if defined(ION_PLATFORM_WINDOWS)
    std::replace(test_case_output.begin(), test_case_output.end(), '/', '\\');
#endif
    EXPECT_EQ(test_case_output, clean_result) << test_case.input;
  }
}

TEST(FileSystem, SplitPath) {
  struct SplitTestCase {
    std::string input;
    std::string first_out;
    std::string second_out;
  };

  const SplitTestCase kSplitTests[] = {
    {{""}, {""}, {""}},       // Empty -> empty
    {{"/"}, {"/"}, {""}},     // Edge case, single slash preserved
    {{"a"}, {""}, {"a"}},     // Edge case, no slash
    {{"/a"}, {"/"}, {"a"}},   // Root slash preserved
    {{"a/b"}, {"a"}, {"b"}},  // Slash between paths consumed in general
    {{"root/folder/file"}, {"root/folder"}, {"file"}},  // Test longer paths
#if defined(ION_PLATFORM_WINDOWS)
    {{"root\\folder\\file"}, {"root\\folder"}, {"file"}},  // Test Windows
#endif
  };

  for (auto& test_case : kSplitTests) {
    auto split_result = FileSystem::SplitPath(test_case.input);
    EXPECT_EQ(test_case.first_out, split_result.first) << test_case.input;
    EXPECT_EQ(test_case.second_out, split_result.second) << test_case.input;
  }
}

TEST(FileSystem, JoinPath) {
  struct JoinPathTestCase {
    JoinPathTestCase(std::string left, std::string right, std::string out)
        : left_in(std::move(left)),
          right_in(std::move(right)),
          output(std::move(out)) {
#if defined(ION_PLATFORM_WINDOWS)
      std::replace(left_in.begin(), left_in.end(), '/', '\\');
      std::replace(right_in.begin(), right_in.end(), '/', '\\');
      std::replace(output.begin(), output.end(), '/', '\\');
#endif
    }
    std::string left_in;
    std::string right_in;
    std::string output;
  };

  const JoinPathTestCase kJoinPathTests[] = {
      {{""}, {""}, {""}},             // Empty, empty -> empty
      {{""}, {"a"}, {"a"}},           // Empty base
      {{"a"}, {""}, {"a/"}},          // Empty rhs
      {{"a"}, {"b"}, {"a/b"}},        // Two folders
      {{"a/"}, {"b"}, {"a/b"}},       // Two folders, first has trailing '/'
      {{"../b"}, {"a"}, {"../b/a"}},  // Multiple '.'
      {{"b//"}, {"c"}, {"b//c"}},     // Multiple '/'
      {{"c"}, {"/a"}, {"/a"}},        // Absolute path: lhs ignored!
  };

  for (auto& test_case : kJoinPathTests) {
    auto join_result =
        FileSystem::JoinPath(test_case.left_in, test_case.right_in);
    // Test the two-arg version.
    EXPECT_EQ(test_case.output, join_result)
        << "L: \"" << test_case.left_in << "\" R: \"" << test_case.right_in
        << "\"";
    // Test the multi-arg version.
    EXPECT_EQ(FileSystem::JoinPath(test_case.output, test_case.output),
              FileSystem::JoinPath(test_case.left_in, test_case.right_in,
                                   test_case.left_in, test_case.right_in))
        << "L: \"" << test_case.left_in << "\" R: \"" << test_case.right_in
        << "\"";
  }
// Test the three-arg version, for kicks.
#if defined(ION_PLATFORM_WINDOWS)
  EXPECT_EQ("foo\\bar\\baz", FileSystem::JoinPath("foo", "bar", "baz"));
#else
  EXPECT_EQ("foo/bar/baz", FileSystem::JoinPath("foo", "bar", "baz"));
#endif
}

}  // namespace
}  // namespace base
}  // namespace seurat
