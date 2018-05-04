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

TEST(JsonValidator, NullOjbectsAreInvalid) {
  Json::Value null_value(Json::nullValue);
  EXPECT_FALSE(JsonValidator::ValidateCapture(null_value).ok());
  EXPECT_FALSE(JsonValidator::ValidateViewGroup(null_value).ok());
  EXPECT_FALSE(JsonValidator::ValidateView(null_value).ok());
  EXPECT_FALSE(JsonValidator::ValidateProjectiveCamera(null_value).ok());
  EXPECT_FALSE(JsonValidator::ValidateDepthImageFile(null_value).ok());
  EXPECT_FALSE(JsonValidator::ValidateImage4File(null_value).ok());
  EXPECT_FALSE(JsonValidator::ValidateImage1File(null_value).ok());
  EXPECT_FALSE(JsonValidator::ValidateMatrix4f(null_value).ok());
}

TEST(JsonValidator, ValidJsonIsValid) {
  std::shared_ptr<base::FileSystem> file_system =
      std::make_shared<base::FileSystem>(testing::GetTestSrcdir());
  const char kJsonFilename[] =
      "com_google_seurat/seurat/ingest/testdata/manifest_v3.json";
  Json::Value root;
  EXPECT_TRUE(
      JsonUtils::ReadJsonFromFile(kJsonFilename, file_system, &root).ok());
  EXPECT_TRUE(JsonValidator::ValidateCapture(root).ok());
}

}  // namespace ingest
}  // namespace seurat
