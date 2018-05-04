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

#ifndef VR_SEURAT_INGEST_JSON_CONVERTER_H_
#define VR_SEURAT_INGEST_JSON_CONVERTER_H_

#include "json/json.h"
#include "seurat/api/api.pb.h"
#include "seurat/ingest/json_utils.h"

namespace seurat {
namespace ingest {

class JsonConverter {
 public:
  // Converts a capture from JSON to proto.
  static void ConvertCapture(const Json::Value& json_capture,
                             api::proto::Capture* proto_capture);

  // Converts a view group from JSON to proto.
  static void ConvertViewGroup(const Json::Value& json_view_group,
                               api::proto::ViewGroup* proto_view_group);

  // Converts a view from JSON to proto.
  static void ConvertView(const Json::Value& json_view,
                          api::proto::View* proto_view);

  // Converts a projective camera from JSON to proto.
  static void ConvertProjectiveCamera(
      const Json::Value& json_projective_camera,
      api::proto::ProjectiveCamera* proto_projective_camera);

  // Converts a depth image file from JSON to proto.
  static void ConvertDepthImageFile(
      const Json::Value& json_depth_image_file,
      api::proto::DepthImageFile* proto_depth_image_file);

  // Converts an image4 file from JSON to proto.
  static void ConvertImage4File(const Json::Value& json_image4_file,
                                api::proto::Image4File* proto_image4_file);

  // Converts an image1 file from JSON to proto.
  static void ConvertImage1File(const Json::Value& json_image1_file,
                                api::proto::Image1File* proto_image1_file);
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_JSON_CONVERTER_H_
