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

#ifndef VR_SEURAT_INGEST_JSON_UTILS_H_
#define VR_SEURAT_INGEST_JSON_UTILS_H_

#include "absl/strings/string_view.h"
#include "json/json.h"
#include "seurat/base/file_system.h"
#include "seurat/base/projective_camera.h"
#include "seurat/base/status.h"

namespace seurat {
namespace ingest {

// This class contains static helper methods for handling Seurat JSON files.
class JsonUtils {
 public:
  // Parses JSON from |json_string| into |root|.
  static base::Status ParseJsonFromString(const std::string& json_string,
                                          Json::Value* root);

  // Reads and parses JSON from |json_file_path| into |root|.
  static base::Status ReadJsonFromFile(
      absl::string_view json_file_path,
      const std::shared_ptr<base::FileSystem>& file_system, Json::Value* root);
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_JSON_UTILS_H_
