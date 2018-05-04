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

#include "seurat/ingest/view_group_loader_factory.h"

#include "gtest/gtest.h"
#include "seurat/base/file_system.h"
#include "seurat/geometry/cube_face.h"
#include "seurat/image/ldi.h"
#include "seurat/testing/test_flags.h"

namespace seurat {
namespace ingest {
namespace {

using base::FileSystem;
using geometry::CubeFace;

TEST(ViewGroupLoaderFactoryTest, CreateFromJson) {
  std::shared_ptr<FileSystem> file_system =
      std::make_shared<FileSystem>(testing::GetTestSrcdir());
  const char kJsonFilename[] =
      "com_google_seurat/seurat/ingest/testdata/manifest_v3.json";
  const int kThreadCount = 1;
  std::unique_ptr<ViewGroupLoader> view_group_loader;
  ViewGroupLoaderFactory::Parameters parameters;
  parameters.thread_count = kThreadCount;
  ViewGroupLoaderFactory factory(parameters, file_system);
  EXPECT_TRUE(
      factory.CreateViewGroupLoader(kJsonFilename, &view_group_loader).ok());
  EXPECT_TRUE(view_group_loader != nullptr);
  EXPECT_EQ(1, view_group_loader->GetNumViewGroups());

  std::vector<std::shared_ptr<base::Camera>> cameras;
  std::vector<image::Ldi4f> ldis;
  EXPECT_TRUE(view_group_loader->LoadViewGroup(0, &cameras, &ldis).ok());
  EXPECT_EQ(1, ldis.size());
  EXPECT_EQ(cameras.size(), ldis.size());
}

TEST(ViewGroupLoaderFactoryTest, CreateFromJson_RootDirectory) {
  // Use a FileSystem with no root directory.
  //
  // If the manifest is specified with an absolute path, then the
  // ViewGroupLoaderFactory should identify the directory containing the
  // manifest to load the LDI from there.
  std::string src_dir = testing::GetTestSrcdir();

  std::shared_ptr<FileSystem> file_system = std::make_shared<FileSystem>("");

  const char kJsonFilename[] =
      "com_google_seurat/seurat/ingest/testdata/manifest_v3.json";
  std::string absolute_manifest_path =
      file_system->JoinPath(src_dir, kJsonFilename);
  const int kThreadCount = 1;
  std::unique_ptr<ViewGroupLoader> view_group_loader;
  ViewGroupLoaderFactory::Parameters parameters;
  parameters.thread_count = kThreadCount;
  ViewGroupLoaderFactory factory(parameters, file_system);
  base::Status status =
      factory.CreateViewGroupLoader(absolute_manifest_path, &view_group_loader);
  EXPECT_TRUE(status.ok());
  EXPECT_TRUE(view_group_loader != nullptr);
  EXPECT_EQ(1, view_group_loader->GetNumViewGroups());

  std::vector<std::shared_ptr<base::Camera>> cameras;
  std::vector<image::Ldi4f> ldis;
  EXPECT_TRUE(view_group_loader->LoadViewGroup(0, &cameras, &ldis).ok());
  EXPECT_EQ(1, ldis.size());
  EXPECT_EQ(cameras.size(), ldis.size());
}

}  // namespace
}  // namespace ingest
}  // namespace seurat
