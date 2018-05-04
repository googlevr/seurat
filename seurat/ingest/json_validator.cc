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

#include "absl/strings/str_cat.h"
#include "seurat/base/projective_camera_util.h"
#include "seurat/geometry/cube_face.h"

namespace seurat {
namespace ingest {

using geometry::CubeFace;

namespace {

struct TypedMember {
  TypedMember(const std::string& name, Json::ValueType type)
      : name(name), type(type) {}
  std::string name;
  Json::ValueType type;
};

base::Status ValidateTypedMember(const TypedMember& typed_member,
                                 const Json::Value& value) {
  Json::ValueType type = typed_member.type;
  std::string name = typed_member.name;
  switch (type) {
    case Json::nullValue:
      if (!value.isNull()) {
        return base::InvalidArgumentError(
            absl::StrCat(name, " is not null: ", value.toStyledString()));
      }
      break;
    case Json::intValue:
      if (!value.isInt()) {
        return base::InvalidArgumentError(
            absl::StrCat(name, " is not an integer: ", value.toStyledString()));
      }
      break;
    case Json::uintValue:
      if (!value.isUInt()) {
        return base::InvalidArgumentError(absl::StrCat(
            name, " is not an unsigned integer: ", value.toStyledString()));
      }
      break;
    case Json::realValue:
      if (!value.isDouble()) {
        return base::InvalidArgumentError(
            absl::StrCat(name, " is not a double: ", value.toStyledString()));
      }
      break;
    case Json::stringValue:
      if (!value.isString()) {
        return base::InvalidArgumentError(
            absl::StrCat(name, " is not a string: ", value.toStyledString()));
      }
      break;
    case Json::booleanValue:
      if (!value.isBool()) {
        return base::InvalidArgumentError(
            absl::StrCat(name, " is not a bool: ", value.toStyledString()));
      }
      break;
    case Json::arrayValue:
      if (!value.isArray()) {
        return base::InvalidArgumentError(
            absl::StrCat(name, " is not an array: ", value.toStyledString()));
      }
      break;
    case Json::objectValue:
      if (!value.isObject()) {
        return base::InvalidArgumentError(
            absl::StrCat(name, " is not an object: ", value.toStyledString()));
      }
      break;
  }

  return base::OkStatus();
}

base::Status ValidateTypedObject(const std::vector<TypedMember>& members,
                                 const Json::Value& value) {
  if (!value.isObject()) {
    return base::InvalidArgumentError(
        absl::StrCat("value is not an object: ", value.toStyledString()));
  }
  for (int i = 0; i < members.size(); ++i) {
    SEURAT_RETURN_IF_ERROR(
        ValidateTypedMember(members[i], value[members[i].name]));
  }
  return base::OkStatus();
}

}  // namespace

base::Status JsonValidator::ValidateCapture(const Json::Value& capture) {
  const TypedMember kViewGroups("view_groups", Json::arrayValue);
  std::vector<TypedMember> members{{kViewGroups}};
  SEURAT_RETURN_IF_ERROR(ValidateTypedObject(members, capture));

  const TypedMember kHeadboxCenter("headbox_center", Json::arrayValue);
  const Json::Value& headbox_center = capture[kHeadboxCenter.name];
  if (!headbox_center.isNull()) {
    // Headbox center is optional. Only validate if it is not null.
    SEURAT_RETURN_IF_ERROR(JsonValidator::ValidatePoint3f(headbox_center));
  }

  const Json::Value& view_groups = capture[kViewGroups.name];
  for (int i = 0; i < view_groups.size(); ++i) {
    SEURAT_RETURN_IF_ERROR(JsonValidator::ValidateViewGroup(view_groups[i]));
  }
  return base::OkStatus();
}

base::Status JsonValidator::ValidateViewGroup(const Json::Value& view_group) {
  const TypedMember kViews("views", Json::arrayValue);
  std::vector<TypedMember> members{{kViews}};
  SEURAT_RETURN_IF_ERROR(ValidateTypedObject(members, view_group));

  const Json::Value& views = view_group[kViews.name];
  for (int i = 0; i < views.size(); ++i) {
    SEURAT_RETURN_IF_ERROR(JsonValidator::ValidateView(views[i]));
  }
  return base::OkStatus();
}

base::Status JsonValidator::ValidateView(const Json::Value& view) {
  const TypedMember kProjectiveCamera("projective_camera", Json::objectValue);
  const TypedMember kDepthImageFile("depth_image_file", Json::objectValue);
  std::vector<TypedMember> members{{kProjectiveCamera, kDepthImageFile}};
  SEURAT_RETURN_IF_ERROR(ValidateTypedObject(members, view));

  SEURAT_RETURN_IF_ERROR(
      JsonValidator::ValidateProjectiveCamera(view[kProjectiveCamera.name]));
  SEURAT_RETURN_IF_ERROR(
      JsonValidator::ValidateDepthImageFile(view[kDepthImageFile.name]));

  return base::OkStatus();
}

base::Status JsonValidator::ValidateProjectiveCamera(
    const Json::Value& camera) {
  const TypedMember kImageWidth("image_width", Json::intValue);
  const TypedMember kImageHeight("image_height", Json::intValue);
  const TypedMember kClipFromEyeMatrix("clip_from_eye_matrix",
                                       Json::arrayValue);
  const TypedMember kWorldFromEyeMatrix("world_from_eye_matrix",
                                        Json::arrayValue);
  const TypedMember kDepthType("depth_type", Json::stringValue);
  std::vector<TypedMember> members{{kImageWidth, kImageHeight,
                                    kClipFromEyeMatrix, kWorldFromEyeMatrix,
                                    kDepthType}};
  SEURAT_RETURN_IF_ERROR(ValidateTypedObject(members, camera));

  SEURAT_RETURN_IF_ERROR(ValidateMatrix4f(camera[kClipFromEyeMatrix.name]));
  SEURAT_RETURN_IF_ERROR(ValidateMatrix4f(camera[kWorldFromEyeMatrix.name]));
  const Json::Value& depth_type = camera[kDepthType.name];
  if (base::DepthTypeFromString(depth_type.asString()) ==
      base::ProjectiveCamera::DepthType::kUnspecified) {
    return base::InvalidArgumentError(
        absl::StrCat("Unknown depth_type: ", depth_type.asString()));
  }

  return base::OkStatus();
}

base::Status JsonValidator::ValidateDepthImageFile(
    const Json::Value& depth_image) {
  std::vector<TypedMember> members{{TypedMember("color", Json::objectValue),
                                    TypedMember("depth", Json::objectValue)}};
  SEURAT_RETURN_IF_ERROR(ValidateTypedObject(members, depth_image));
  return base::OkStatus();
}

base::Status JsonValidator::ValidateImage4File(const Json::Value& image4_file) {
  std::vector<TypedMember> members{
      {TypedMember("path", Json::stringValue),
       TypedMember("channel_0", Json::stringValue),
       TypedMember("channel_1", Json::stringValue),
       TypedMember("channel_2", Json::stringValue),
       TypedMember("channel_alpha", Json::stringValue)}};
  SEURAT_RETURN_IF_ERROR(ValidateTypedObject(members, image4_file));
  return base::OkStatus();
}

base::Status JsonValidator::ValidateImage1File(const Json::Value& image1_file) {
  std::vector<TypedMember> members{
      {TypedMember("path", Json::stringValue),
       TypedMember("channel_0", Json::stringValue)}};
  SEURAT_RETURN_IF_ERROR(ValidateTypedObject(members, image1_file));
  return base::OkStatus();
}

base::Status JsonValidator::ValidateMatrix4f(const Json::Value& matrix) {
  if (!matrix.isArray()) {
    return base::InvalidArgumentError(
        absl::StrCat("matrix is not an array: ", matrix.toStyledString()));
  }

  if (matrix.size() != 16) {
    return base::InvalidArgumentError(
        absl::StrCat("size of matrix != 16: ", matrix.toStyledString()));
  }

  for (int i = 0; i < 16; ++i) {
    if (!matrix[i].isNumeric()) {
      return base::InvalidArgumentError(absl::StrCat(
          "matrix[", i, "] is not numeric: ", matrix[i].toStyledString()));
    }
  }

  return base::OkStatus();
}

base::Status JsonValidator::ValidatePoint3f(const Json::Value& point) {
  if (!point.isArray()) {
    return base::InvalidArgumentError(
        absl::StrCat("point is not an array: ", point.toStyledString()));
  }

  if (point.size() != 3) {
    return base::InvalidArgumentError(
        absl::StrCat("size of point != 3: ", point.toStyledString()));
  }

  for (int i = 0; i < 3; ++i) {
    if (!point[i].isNumeric()) {
      return base::InvalidArgumentError(absl::StrCat(
          "point[", i, "] is not numeric: ", point[i].toStyledString()));
    }
  }

  return base::OkStatus();
}

}  // namespace ingest
}  // namespace seurat
