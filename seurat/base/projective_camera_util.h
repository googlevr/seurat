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

#ifndef VR_SEURAT_BASE_PROJECTIVE_CAMERA_UTIL_H_
#define VR_SEURAT_BASE_PROJECTIVE_CAMERA_UTIL_H_

#include "absl/strings/string_view.h"
#include "json/json.h"
#include "seurat/api/camera.pb.h"
#include "seurat/base/projective_camera.h"
#include "seurat/base/status.h"

namespace seurat {
namespace base {

// Returns the depth type for the given string.
base::ProjectiveCamera::DepthType DepthTypeFromString(
    absl::string_view depth_type_string);

// Returns a string representation of the |depth_type|.
std::string StringFromDepthType(base::ProjectiveCamera::DepthType depth_type);

// Returns a ProjectiveCamera read from protocol buffer object |proto|.
ProjectiveCamera ProjectiveCameraFromProto(
    const api::proto::ProjectiveCamera& proto);

// Writes |camera| to protocol buffer object |proto|.
void ProjectiveCameraToProto(const ProjectiveCamera& camera,
                             api::proto::ProjectiveCamera* proto);

// Returns a JSON serialization of the |camera|.
Json::Value JsonFromProjectiveCamera(const ProjectiveCamera& camera);

// Deserialized |projective_camera| into a ProjectiveCamera and returns it in a
// shared pointer. Assumes that |projective_camera| is a valid JSON
// serialization of a projective camera.
base::ProjectiveCamera ProjectiveCameraFromJson(
    const Json::Value& projective_camera);

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_PROJECTIVE_CAMERA_UTIL_H_
