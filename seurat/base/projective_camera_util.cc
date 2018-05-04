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

#include "seurat/base/projective_camera_util.h"

#include "ion/math/matrix.h"
#include "ion/math/matrixutils.h"
#include "ion/math/vector.h"
#include "absl/strings/str_cat.h"
#include "seurat/api/camera.pb.h"
#include "seurat/base/math_util.h"

namespace seurat {
namespace base {

using ion::math::Matrix4f;
using ion::math::Vector2i;

ProjectiveCamera::DepthType DepthTypeFromString(
    absl::string_view depth_type_string) {
  ProjectiveCamera::DepthType depth_type =
      ProjectiveCamera::DepthType::kUnspecified;
  // TODO(ernstm): Define all string constants in a central location.
  if (depth_type_string == "WINDOW_Z") {
    depth_type = ProjectiveCamera::DepthType::kWindowZ;
  } else if (depth_type_string == "EYE_Z") {
    depth_type = ProjectiveCamera::DepthType::kEyeZ;
  } else if (depth_type_string == "RAY_DEPTH") {
    depth_type = ProjectiveCamera::DepthType::kRayDepth;
  } else {
    // Log an error (and return kUnspecified), if the depth type is unknown.
    LOG(ERROR) << "Unknown depth type:" << depth_type_string;
  }
  return depth_type;
}

std::string StringFromDepthType(ProjectiveCamera::DepthType depth_type) {
  switch (depth_type) {
    case ProjectiveCamera::DepthType::kUnspecified:
      return "UNSPECIFIED";
    case ProjectiveCamera::DepthType::kWindowZ:
      return "WINDOW_Z";
    case ProjectiveCamera::DepthType::kEyeZ:
      return "EYE_Z";
    case ProjectiveCamera::DepthType::kRayDepth:
      return "RAY_DEPTH";
  }
}

ProjectiveCamera ProjectiveCameraFromProto(
    const api::proto::ProjectiveCamera& proto) {
  const ion::math::Vector2i image_size(proto.image_size().x(),
                                       proto.image_size().y());
  const ion::math::Matrix4f world_from_eye =
      Matrix4fFromProto(proto.world_from_eye());
  const ion::math::Matrix4f eye_from_world = ion::math::Inverse(world_from_eye);
  const ion::math::Matrix4f clip_from_eye(
      Matrix4fFromProto(proto.clip_from_eye()));
  // TODO(davejr, b/36684150) - Reduce the projective camera constructor back to
  // one parameter - just the camera matrix (and fix all the test code this
  // requires).
  const ProjectiveCamera::DepthType depth_type =
      static_cast<ProjectiveCamera::DepthType>(proto.depth_type());
  return ProjectiveCamera(image_size, clip_from_eye, eye_from_world,
                          world_from_eye, depth_type);
}

void ProjectiveCameraToProto(const ProjectiveCamera& camera,
                             api::proto::ProjectiveCamera* proto) {
  proto->mutable_image_size()->set_x(camera.GetImageSize()[0]);
  proto->mutable_image_size()->set_y(camera.GetImageSize()[1]);
  Matrix4fToProto(camera.GetWorldFromEye(), proto->mutable_world_from_eye());
  Matrix4fToProto(camera.GetClipFromEye(), proto->mutable_clip_from_eye());
  proto->set_depth_type(
      static_cast<api::proto::DepthType>(camera.GetDepthType()));
}

Json::Value JsonFromProjectiveCamera(const ProjectiveCamera& camera) {
  Json::Value value(Json::objectValue);
  value["image_width"] = camera.GetImageSize()[0];
  value["image_height"] = camera.GetImageSize()[1];
  value["clip_from_eye_matrix"] = JsonFromMatrix4f(camera.GetClipFromEye());
  value["world_from_eye_matrix"] = JsonFromMatrix4f(camera.GetWorldFromEye());
  value["depth_type"] = StringFromDepthType(camera.GetDepthType());
  return value;
}

ProjectiveCamera ProjectiveCameraFromJson(const Json::Value& camera) {
  const Json::Value& width = camera["image_width"];
  const Json::Value& height = camera["image_height"];
  const Vector2i image_size(width.asInt(), height.asInt());
  const Matrix4f clip_from_eye =
      Matrix4fFromJson(camera["clip_from_eye_matrix"]);
  const Matrix4f eye_from_world =
      ion::math::Inverse(Matrix4fFromJson(camera["world_from_eye_matrix"]));
  const ProjectiveCamera::DepthType depth_type =
      DepthTypeFromString(camera["depth_type"].asString());
  // TODO(davejr, b/36684150) - Reduce the projective camera constructor back to
  // one parameter - just the camera matrix (and fix all the test code this
  // requires).
  return ProjectiveCamera(image_size, clip_from_eye, eye_from_world,
                          depth_type);
}

}  // namespace base
}  // namespace seurat
