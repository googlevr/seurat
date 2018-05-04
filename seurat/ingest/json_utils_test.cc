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

#include "seurat/ingest/json_validator.h"

#include "gtest/gtest.h"
#include "json/json.h"
#include "seurat/base/file_system.h"
#include "seurat/ingest/json_utils.h"
#include "seurat/testing/test_flags.h"

namespace seurat {
namespace ingest {

TEST(JsonUtils, ReadingValidJsonSucceeds) {
  std::shared_ptr<base::FileSystem> file_system =
      std::make_shared<base::FileSystem>(testing::GetTestSrcdir());
  const char kJsonFilename[] =
      "com_google_seurat/seurat/ingest/testdata/manifest_v3.json";
  Json::Value root;
  EXPECT_TRUE(
      JsonUtils::ReadJsonFromFile(kJsonFilename, file_system, &root).ok());
}

TEST(JsonUtils, ReadingInvalidJsonFails) {
  std::shared_ptr<base::FileSystem> file_system =
      std::make_shared<base::FileSystem>(testing::GetTestSrcdir());
  const char kJsonFilename[] =
      "com_google_seurat/seurat/ingest/testdata/invalid_manifest.json";
  Json::Value root;
  EXPECT_FALSE(
      JsonUtils::ReadJsonFromFile(kJsonFilename, file_system, &root).ok());
  EXPECT_TRUE(root.isNull());
}

}  // namespace ingest
}  // namespace seurat
