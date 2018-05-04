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

#include "seurat/component/ice_io.h"

#include <memory>

#include "gtest/gtest.h"
#include "absl/strings/string_view.h"
#include "seurat/base/file_system.h"
#include "seurat/base/structured_io.h"
#include "seurat/component/component_tests_utils.h"
#include "seurat/testing/test_flags.h"

namespace seurat {
namespace component {
namespace {

using base::ArrayByteSource;
using base::FileSystem;
using base::StringByteSink;
using base::StructureSink;
using base::StructureSource;

base::Status WriteIceFile(const Component& component, FileSystem* file_system,
                          absl::string_view filename) {
  std::string buf;
  StringByteSink ice_string_byte_sink(&buf);
  StructureSink ice_sink(&ice_string_byte_sink);
  WriteIce(component, &ice_sink);
  return file_system->SetContents(filename, buf);
}

base::Status ReadIceFile(absl::string_view filename, FileSystem* file_system,
                         std::unique_ptr<const Component>* component) {
  // Read the component hierarchy from the file.
  std::string file_string;
  base::Status status = file_system->GetContents(filename, &file_string);
  if (status.ok()) {
    ArrayByteSource byte_source(file_string);
    StructureSource ice_source(&byte_source);
    *component = ReadIce(&ice_source);
  }

  return status;
}

// Reads and writes a simple nested hierarchy of Component instances.
TEST(IceIoTest, WriteReadIce) {
  FakeComponent source_component(
      std::unique_ptr<FakeRenderable>(new FakeRenderable));

  base::FileSystem file_system(testing::GetTestTmpdir());
  const std::string ice_filename("test_ice_header.ice");

  base::Status write_status =
      WriteIceFile(source_component, &file_system, ice_filename);
  ASSERT_TRUE(write_status.ok())
      << "Writing output \"" << ice_filename << "\" failed";

  std::unique_ptr<const Component> loaded_ice_component;
  base::Status read_status =
      ReadIceFile(ice_filename, &file_system, &loaded_ice_component);
  ASSERT_TRUE(read_status.ok())
      << "Reading generated ice file \"" << ice_filename << "\" failed";
  ASSERT_NE(loaded_ice_component.get(), nullptr);

  EXPECT_EQ(source_component, *loaded_ice_component);
}

// Builds an ice file of the prior API revision, and verifies reads using the
// current API fail.
TEST(IceIoTest, TestPriorIceVersion) {
  std::string buf;
  StringByteSink ice_file_bytes(&buf);
  StructureSink ice_sink(&ice_file_bytes);
  FakeComponent source_component(
      std::unique_ptr<FakeRenderable>(new FakeRenderable));

  static const int kOlderFormat = ice_io_internal::kCurrent - 1;
  ice_io_internal::IceFormat test_ice_format(kOlderFormat);
  test_ice_format.WriteIce(source_component, &ice_sink);

  // Get the bytes.
  size_t ice_file_buffer_size = buf.size();
  ASSERT_GT(ice_file_buffer_size, 0);
  absl::string_view ice_bytes(buf.data(), ice_file_buffer_size);

  // Attempt to read the data with the current API.
  ArrayByteSource file_stream(ice_bytes);
  StructureSource ice_source(&file_stream);
  auto loaded_ice_component = ReadIce(&ice_source);
  EXPECT_EQ(loaded_ice_component.get(), nullptr);
}

// Builds an ICE file with a newer API revision, and verifies reads using the
// current API fail.
TEST(IceIoTest, TestNextIceVersion) {
  // Build an ice file of a newer revision than current.
  std::string buf;
  StringByteSink strsink(&buf);
  StructureSink ice_sink(&strsink);
  FakeComponent source_component(
      std::unique_ptr<FakeRenderable>(new FakeRenderable));

  static const int kNewerFormat = ice_io_internal::kCurrent + 1;
  ice_io_internal::IceFormat test_ice_format(kNewerFormat);
  test_ice_format.WriteIce(source_component, &ice_sink);

  // Get the bytes.
  size_t ice_file_buffer_size = buf.size();
  ASSERT_GT(ice_file_buffer_size, 0);
  absl::string_view ice_bytes(buf.data(), ice_file_buffer_size);

  // Attempt to read the data with an older format version.
  ArrayByteSource file_stream(ice_bytes);
  StructureSource ice_source(&file_stream);
  auto loaded_ice_component = ReadIce(&ice_source);
  EXPECT_EQ(loaded_ice_component.get(), nullptr);
}

// Tests deserializing various broken (truncated) ICE format byte sequences
// fails.
//
// NOTE! This test is disabled because the Structure IO code currently CHECKs on
// all read errors.
TEST(IceIoTest, DISABLED_TestTruncatedIce) {
  // Build an ice file.
  std::string buf;
  StringByteSink strsink(&buf);
  StructureSink ice_sink(&strsink);
  FakeComponent source_component(
      std::unique_ptr<FakeRenderable>(new FakeRenderable));

  WriteIce(source_component, &ice_sink);

  // Get the bytes.
  size_t ice_file_size = buf.size();
  ASSERT_GT(ice_file_size, 0);

  // Test truncating the stream at each point results in a failure to load a
  // component.
  for (int test_truncation_size = 0; test_truncation_size <= ice_file_size;
       ++test_truncation_size) {
    std::cout << "Testing trunc @ " << test_truncation_size << std::endl;
    ArrayByteSource truncated_file_stream(
        absl::string_view(buf.data(), test_truncation_size));
    StructureSource ice_source(&truncated_file_stream);
    auto loaded_ice_component = ReadIce(&ice_source);
    if (test_truncation_size < ice_file_size) {
      EXPECT_EQ(loaded_ice_component.get(), nullptr);
    } else {
      // The untruncated read should succeed!
      EXPECT_NE(loaded_ice_component.get(), nullptr);
    }
  }
}

}  // namespace
}  // namespace component
}  // namespace seurat
