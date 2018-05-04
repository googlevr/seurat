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

#include <array>
#include <string>

#include "ion/math/matrix.h"
#include "ion/math/transformutils.h"
#include "gtest/gtest.h"
#include "seurat/base/projective_camera.h"

namespace seurat {
namespace base {
namespace {

using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2i;
using ion::math::Vector3f;

TEST(ProjectiveCamera, DepthTypesConvertCorrectly) {
  const ProjectiveCamera::DepthType kDepthTypes[] = {
      ProjectiveCamera::DepthType::kUnspecified,
      ProjectiveCamera::DepthType::kWindowZ, ProjectiveCamera::DepthType::kEyeZ,
      ProjectiveCamera::DepthType::kRayDepth};
  for (const auto depth_type : kDepthTypes) {
    // The switch causes the compiler to issue an error if a newly added enum
    // value is not handled by this test.
    switch (depth_type) {
      case ProjectiveCamera::DepthType::kUnspecified:
        EXPECT_EQ("UNSPECIFIED", StringFromDepthType(depth_type));
        EXPECT_EQ(depth_type, DepthTypeFromString("UNSPECIFIED"));
        EXPECT_EQ(static_cast<int>(depth_type),
                  static_cast<int>(api::proto::DEPTH_TYPE_UNSPECIFIED));
        break;
      case ProjectiveCamera::DepthType::kWindowZ:
        EXPECT_EQ("WINDOW_Z", StringFromDepthType(depth_type));
        EXPECT_EQ(depth_type, DepthTypeFromString("WINDOW_Z"));
        EXPECT_EQ(static_cast<int>(depth_type),
                  static_cast<int>(api::proto::DEPTH_TYPE_WINDOW_Z));
        break;
      case ProjectiveCamera::DepthType::kEyeZ:
        EXPECT_EQ("EYE_Z", StringFromDepthType(depth_type));
        EXPECT_EQ(depth_type, DepthTypeFromString("EYE_Z"));
        EXPECT_EQ(static_cast<int>(depth_type),
                  static_cast<int>(api::proto::DEPTH_TYPE_EYE_Z));
        break;
      case ProjectiveCamera::DepthType::kRayDepth:
        EXPECT_EQ("RAY_DEPTH", StringFromDepthType(depth_type));
        EXPECT_EQ(depth_type, DepthTypeFromString("RAY_DEPTH"));
        EXPECT_EQ(static_cast<int>(depth_type),
                  static_cast<int>(api::proto::DEPTH_TYPE_RAY_DEPTH));
        break;
    }
  }
}

// Verify serialization and deserialization based on protos.
TEST(ProjectiveCameraUtil, SerializeAndDeserializeProto) {
  const float kNear = 0.1f;
  const float kFar = 100.0f;
  const Vector2i kImageSize(32, 24);
  const std::array<ProjectiveCamera::DepthType, 4> kDepthTypes{
      {ProjectiveCamera::DepthType::kUnspecified,
       ProjectiveCamera::DepthType::kWindowZ,
       ProjectiveCamera::DepthType::kEyeZ,
       ProjectiveCamera::DepthType::kRayDepth}};
  for (const auto& depth_type : kDepthTypes) {
    const ProjectiveCamera camera(
        kImageSize,
        ion::math::PerspectiveMatrixFromFrustum(-kNear, kNear, -kNear, kNear,
                                                kNear, kFar),
        ion::math::Matrix4f::Identity(), depth_type);
    api::proto::ProjectiveCamera proto;
    ProjectiveCameraToProto(camera, &proto);
    std::string proto_string = proto.SerializeAsString();
    api::proto::ProjectiveCamera deserialized_proto;
    deserialized_proto.ParseFromString(proto_string);
    ProjectiveCamera deserialized =
        ProjectiveCameraFromProto(deserialized_proto);
    EXPECT_EQ(camera, deserialized);
  }
}

TEST(ProjectiveCameraUtil, SerializeAndDeserializeJson) {
  const float kNear = 0.1f;
  const float kFar = 100.0f;
  const Vector2i kImageSize(32, 24);
  const std::array<ProjectiveCamera::DepthType, 4> kDepthTypes{
      {ProjectiveCamera::DepthType::kUnspecified,
       ProjectiveCamera::DepthType::kWindowZ,
       ProjectiveCamera::DepthType::kEyeZ,
       ProjectiveCamera::DepthType::kRayDepth}};
  for (const auto& depth_type : kDepthTypes) {
    const ProjectiveCamera camera(
        kImageSize,
        ion::math::PerspectiveMatrixFromFrustum(-kNear, kNear, -kNear, kNear,
                                                kNear, kFar),
        ion::math::Matrix4f::Identity(), depth_type);
    Json::Value camera_json = JsonFromProjectiveCamera(camera);
    const ProjectiveCamera camera_copy = ProjectiveCameraFromJson(camera_json);
    EXPECT_EQ(camera, camera_copy);
  }
}
}  // namespace
}  // namespace base
}  // namespace seurat
