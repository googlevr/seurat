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

#include "seurat/ingest/json_utils.h"

#include "seurat/base/projective_camera_util.h"
#include "seurat/ingest/json_validator.h"

namespace seurat {
namespace ingest {

base::Status JsonUtils::ParseJsonFromString(const std::string& json_string,
                                            Json::Value* root) {
  Json::Reader reader;
  if (!reader.parse(json_string, *root)) {
    *root = Json::nullValue;
    return base::FailedPreconditionError(reader.getFormattedErrorMessages());
  }
  return base::OkStatus();
}

base::Status JsonUtils::ReadJsonFromFile(
    absl::string_view json_file_path,
    const std::shared_ptr<base::FileSystem>& file_system, Json::Value* root) {
  std::string json_string;
  SEURAT_RETURN_IF_ERROR(
      file_system->GetContents(json_file_path, &json_string));
  SEURAT_RETURN_IF_ERROR(ParseJsonFromString(json_string, root));
  return base::OkStatus();
}

}  // namespace ingest
}  // namespace seurat
