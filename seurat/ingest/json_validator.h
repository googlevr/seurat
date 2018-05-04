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

#ifndef VR_SEURAT_INGEST_JSON_VALIDATOR_H_
#define VR_SEURAT_INGEST_JSON_VALIDATOR_H_

#include "json/json.h"
#include "seurat/base/status.h"

namespace seurat {
namespace ingest {

class JsonValidator {
 public:
  // Validates that |capture| is a correct JSON serialization of a Seurat
  // capture.
  static base::Status ValidateCapture(const Json::Value& capture);

  // Validates if |view_group| is a correct JSON serialization of a view group.
  static base::Status ValidateViewGroup(const Json::Value& view_group);

  // Validates if |view| is a correct JSON serialization of a view.
  static base::Status ValidateView(const Json::Value& view);

  // Validates if |projective_camera| is a correct JSON serialization of a
  // ProjectiveCamera.
  static base::Status ValidateProjectiveCamera(const Json::Value& camera);

  // Validates if |depth_image| is a correct JSON serialization of a depth
  // image file.
  static base::Status ValidateDepthImageFile(const Json::Value& depth_image);

  // Validates if |image_file| is a correct JSON serialization of a four-channel
  // image file.
  static base::Status ValidateImage4File(const Json::Value& image4_file);

  // Validates if |image_file| is a correct JSON serialization of a one-channel
  // image file.
  static base::Status ValidateImage1File(const Json::Value& image1_file);

  // Validates if |matrix| is a correct JSON serialization of a Matrix4f.
  static base::Status ValidateMatrix4f(const Json::Value& matrix);

  // Validates if |point| is a correct JSON serialization of a Point3f.
  static base::Status ValidatePoint3f(const Json::Value& point);
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_JSON_VALIDATOR_H_
