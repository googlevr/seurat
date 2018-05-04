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

#include "seurat/base/camera_util.h"

#include <array>
#include <string>

#include "ion/math/matrix.h"
#include "ion/math/transformutils.h"
#include "gtest/gtest.h"
#include "seurat/base/projective_camera.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace base {
namespace {

using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2i;
using ion::math::Vector3f;

TEST(CameraUtil, DepthTypesConvertCorrectly) {
  const float kNear = 0.1f;
  const float kFar = 100.0f;
  const Vector2i kImageSize(32, 24);
  const ProjectiveCamera::DepthType kDepthType =
      ProjectiveCamera::DepthType::kEyeZ;

  auto original_camera = std::make_shared<ProjectiveCamera>(
      kImageSize,
      ion::math::PerspectiveMatrixFromFrustum(-kNear, kNear, -kNear, kNear,
                                              kNear, kFar),
      ion::math::Matrix4f::Identity(), kDepthType);

  Vector3f translation(2.0f, 8.0f, -6.0f);
  std::unique_ptr<Camera> translated_camera =
      TranslateCamera(translation, original_camera);

  EXPECT_EQ(original_camera->GetImageSize(), translated_camera->GetImageSize());

  Point3f original_eye = ion::math::ProjectPoint(
      original_camera->GetWorldFromEye(), Point3f::Zero());
  Point3f translated_eye = ion::math::ProjectPoint(
      translated_camera->GetWorldFromEye(), Point3f::Zero());
  EXPECT_VECTOR_NEAR(original_eye + translation, translated_eye, 1e-5f);

  EXPECT_VECTOR_NEAR(original_camera->RayOrigin({8, 3}) + translation,
                     translated_camera->RayOrigin({8, 3}), 1e-5f);

  EXPECT_VECTOR_NEAR(original_camera->RayDirection({3, 7}),
                     translated_camera->RayDirection({3, 7}), 1e-5f);

  EXPECT_VECTOR_NEAR(original_camera->RayEnd({23, 31}, 42.0f) + translation,
                     translated_camera->RayEnd({23, 31}, 42.0f), 1e-5f);
}

}  // namespace
}  // namespace base
}  // namespace seurat
