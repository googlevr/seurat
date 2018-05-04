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

#include "seurat/ingest/json_converter.h"

#include "json/json.h"
#include "seurat/api/camera.pb.h"
#include "seurat/api/image.pb.h"
#include "seurat/api/math.pb.h"
#include "seurat/base/projective_camera_util.h"

namespace seurat {
namespace ingest {

using seurat::base::ProjectiveCamera;

void JsonConverter::ConvertCapture(const Json::Value& json_capture,
                                   api::proto::Capture* proto_capture) {
  const Json::Value json_headbox_center = json_capture["headbox_center"];
  if (!json_headbox_center.isNull()) {
    proto_capture->mutable_headbox_center()->set_x(
        json_headbox_center[0].asFloat());
    proto_capture->mutable_headbox_center()->set_y(
        json_headbox_center[1].asFloat());
    proto_capture->mutable_headbox_center()->set_z(
        json_headbox_center[2].asFloat());
  }

  const Json::Value& json_view_groups = json_capture["view_groups"];
  for (int i = 0; i < json_view_groups.size(); ++i) {
    api::proto::ViewGroup* proto_view_group = proto_capture->add_view_groups();
    ConvertViewGroup(json_view_groups[i], proto_view_group);
  }
}

void JsonConverter::ConvertViewGroup(const Json::Value& json_view_group,
                                     api::proto::ViewGroup* proto_view_group) {
  const Json::Value& json_views = json_view_group["views"];
  for (int i = 0; i < json_views.size(); ++i) {
    api::proto::View* proto_view = proto_view_group->add_views();
    ConvertView(json_views[i], proto_view);
  }
}

void JsonConverter::ConvertView(const Json::Value& json_view,
                                api::proto::View* proto_view) {
  ConvertProjectiveCamera(json_view["projective_camera"],
                          proto_view->mutable_camera()->mutable_projective());
  ConvertDepthImageFile(json_view["depth_image_file"],
                        proto_view->mutable_ldi()->mutable_depth_image_file());
}

void JsonConverter::ConvertProjectiveCamera(
    const Json::Value& json_projective_camera,
    api::proto::ProjectiveCamera* proto_projective_camera) {
  // TODO(ernstm): Do a direct conversion here?
  base::ProjectiveCamera camera =
      base::ProjectiveCameraFromJson(json_projective_camera);
  ProjectiveCameraToProto(camera, proto_projective_camera);
}

void JsonConverter::ConvertDepthImageFile(
    const Json::Value& json_depth_image_file,
    api::proto::DepthImageFile* proto_depth_image_file) {
  ConvertImage4File(json_depth_image_file["color"],
                    proto_depth_image_file->mutable_color());
  ConvertImage1File(json_depth_image_file["depth"],
                    proto_depth_image_file->mutable_depth());
}

void JsonConverter::ConvertImage4File(
    const Json::Value& json_image4_file,
    api::proto::Image4File* proto_image4_file) {
  proto_image4_file->set_path(json_image4_file["path"].asString());
  proto_image4_file->set_channel_0(json_image4_file["channel_0"].asString());
  proto_image4_file->set_channel_1(json_image4_file["channel_1"].asString());
  proto_image4_file->set_channel_2(json_image4_file["channel_2"].asString());
  proto_image4_file->set_channel_alpha(
      json_image4_file["channel_alpha"].asString());
}

void JsonConverter::ConvertImage1File(
    const Json::Value& json_image1_file,
    api::proto::Image1File* proto_image1_file) {
  proto_image1_file->set_path(json_image1_file["path"].asString());
  proto_image1_file->set_channel_0(json_image1_file["channel_0"].asString());
}

}  // namespace ingest
}  // namespace seurat
