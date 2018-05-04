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

#include "seurat/base/projective_camera.h"

#include <random>

#include "ion/math/matrix.h"
#include "ion/math/matrixutils.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "gtest/gtest.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace base {

using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2i;
using ion::math::Vector3f;

const int kWidth = 32;
const int kHeight = 24;
const Vector2i kImageSize = {kWidth, kHeight};
const float kLeft = -0.8f;
const float kRight = 1.0f;
const float kBottom = -0.7f;
const float kTop = 0.9f;
const float kNear = 1.0f;
const float kFar = 100.0f;
const Vector3f kTranslation(3.0f, -1.0f, 2.0f);
const float kFarScale = kFar / kNear;
// clang-format off
// World-space positions of the 8 frustum corners.
const Point3f kLeftBottomNear =
    Point3f(kLeft, kBottom, -kNear) - kTranslation;
const Point3f kLeftBottomFar =
    Point3f(kLeft * kFarScale, kBottom * kFarScale, -kFar) - kTranslation;
const Point3f kLeftTopNear =
    Point3f(kLeft, kTop, -kNear) - kTranslation;
const Point3f kLeftTopFar =
    Point3f(kLeft * kFarScale, kTop * kFarScale, -kFar) - kTranslation;
const Point3f kRightBottomNear =
    Point3f(kRight, kBottom, -kNear) - kTranslation;
const Point3f kRightBottomFar =
    Point3f(kRight * kFarScale, kBottom * kFarScale, -kFar) - kTranslation;
const Point3f kRightTopNear =
    Point3f(kRight, kTop, -kNear) - kTranslation;
const Point3f kRightTopFar =
    Point3f(kRight * kFarScale, kTop * kFarScale, -kFar) - kTranslation;
// clang-format on

ProjectiveCamera CreateTestCamera(ProjectiveCamera::DepthType depth_type =
                                      ProjectiveCamera::DepthType::kWindowZ) {
  const Matrix4f kClipFromEyeMatrix(ion::math::PerspectiveMatrixFromFrustum(
      kLeft, kRight, kBottom, kTop, kNear, kFar));
  const Matrix4f kEyeFromWorldMatrix(
      ion::math::TranslationMatrix(kTranslation));
  return ProjectiveCamera(kImageSize, kClipFromEyeMatrix, kEyeFromWorldMatrix,
                          depth_type);
}

TEST(ProjectiveCamera, Constructor) {
  const Matrix4f clip_from_eye_matrix(ion::math::PerspectiveMatrixFromFrustum(
      kLeft, kRight, kBottom, kTop, kNear, kFar));
  const Matrix4f eye_from_world_matrix(
      ion::math::TranslationMatrix(kTranslation));
  const ProjectiveCamera camera(kImageSize, clip_from_eye_matrix,
                                eye_from_world_matrix);
  EXPECT_EQ(kImageSize, camera.GetImageSize());
  EXPECT_EQ(clip_from_eye_matrix, camera.GetClipFromEye());
  EXPECT_EQ(eye_from_world_matrix, camera.GetEyeFromWorld());

  const float w_2 = kImageSize[0] * 0.5f;
  const float h_2 = kImageSize[1] * 0.5f;
  // clang-format off
  const Matrix4f window_from_clip_matrix(w_2,  0.0f, 0.0f, w_2,
                                         0.0f, h_2,  0.0f, h_2,
                                         0.0f, 0.0f, 0.5f, 0.5f,
                                         0.0f, 0.0f, 0.0f, 1.0f);
  // clang-format on
  const Matrix4f world_from_eye_matrix =
      ion::math::Inverse(eye_from_world_matrix);
  const Matrix4f window_from_eye_matrix =
      window_from_clip_matrix * clip_from_eye_matrix;
  const Matrix4f eye_from_window_matrix =
      ion::math::Inverse(window_from_eye_matrix);
  const Matrix4f window_from_world_matrix =
      window_from_eye_matrix * eye_from_world_matrix;
  const Matrix4f world_from_window_matrix =
      world_from_eye_matrix * eye_from_window_matrix;
  const Matrix4f eye_from_clip_matrix =
      ion::math::Inverse(clip_from_eye_matrix);
  EXPECT_EQ(world_from_eye_matrix, camera.GetWorldFromEye());
  EXPECT_EQ(world_from_eye_matrix, camera.GetWorldFromEye());
  EXPECT_EQ(window_from_eye_matrix, camera.GetWindowFromEye());
  EXPECT_EQ(eye_from_window_matrix, camera.GetEyeFromWindow());
  EXPECT_EQ(window_from_world_matrix, camera.GetWindowFromWorld());
  EXPECT_EQ(world_from_window_matrix, camera.GetWorldFromWindow());
  EXPECT_EQ(eye_from_clip_matrix, camera.GetEyeFromClip());
}

TEST(ProjectiveCamera, EqualCamerasAreEqual) {
  const ProjectiveCamera camera_a = CreateTestCamera();
  const ProjectiveCamera camera_b = CreateTestCamera();
  EXPECT_EQ(camera_a, camera_b);
}

TEST(ProjectiveCamera, DifferentCamerasAreNotEqual) {
  Matrix4f clip_from_eye_matrix(ion::math::PerspectiveMatrixFromFrustum(
      kLeft, kRight, kBottom, kTop, kNear, kFar));
  Matrix4f eye_from_world_matrix(ion::math::TranslationMatrix(kTranslation));
  ProjectiveCamera camera_a(kImageSize, clip_from_eye_matrix,
                            eye_from_world_matrix);
  ProjectiveCamera camera_b(kImageSize + Vector2i{1, 0}, clip_from_eye_matrix,
                            eye_from_world_matrix);
  ProjectiveCamera camera_c(kImageSize + Vector2i{0, 1}, clip_from_eye_matrix,
                            eye_from_world_matrix);
  ProjectiveCamera camera_d(kImageSize,
                            clip_from_eye_matrix + Matrix4f::Identity(),
                            eye_from_world_matrix);
  ProjectiveCamera camera_e(kImageSize, clip_from_eye_matrix,
                            eye_from_world_matrix + Matrix4f::Identity());
  EXPECT_NE(camera_a, camera_b);
  EXPECT_NE(camera_a, camera_d);
  EXPECT_NE(camera_a, camera_d);
  EXPECT_NE(camera_a, camera_e);
}

TEST(ProjectiveCamera, RayOrigin) {
  ProjectiveCamera camera = CreateTestCamera();

  // Verify that the 4 corners of the image plane map to the four near-plane
  // frustum corners.
  EXPECT_VECTOR_NEAR(kLeftBottomNear,
                     camera.RayOriginFloat(Point2f(0.0f, 0.0f)), 1e-5f);
  EXPECT_VECTOR_NEAR(kLeftTopNear,
                     camera.RayOriginFloat(Point2f(0.0f, kHeight)), 1e-5f);
  EXPECT_VECTOR_NEAR(kRightBottomNear,
                     camera.RayOriginFloat(Point2f(kWidth, 0.0f)), 1e-5f);
  EXPECT_VECTOR_NEAR(kRightTopNear,
                     camera.RayOriginFloat(Point2f(kWidth, kHeight)), 1e-5f);

  // Test integer version.
  EXPECT_EQ(camera.RayOriginFloat(Point2f(0.5f, 0.5f)),
            camera.RayOrigin(Point2i(0, 0)));
  EXPECT_EQ(camera.RayOriginFloat(Point2f(kWidth - 0.5f, kHeight - 0.5f)),
            camera.RayOrigin(Point2i(kWidth - 1, kHeight - 1)));
}

TEST(ProjectiveCamera, RayDirection) {
  ProjectiveCamera camera = CreateTestCamera();
  Point3f eye_position_world =
      camera.GetWorldFromEye() * ion::math::Point3f::Zero();

  // Verify that the 8 corner points of the frustum are handled correctly.
  EXPECT_VECTOR_NEAR(
      ion::math::Normalized(kLeftBottomNear - eye_position_world),
      ion::math::Normalized(camera.RayDirectionFloat(Point2f(0.0f, 0.0f))),
      1e-5f);
  EXPECT_VECTOR_NEAR(
      ion::math::Normalized(kLeftTopNear - eye_position_world),
      ion::math::Normalized(camera.RayDirectionFloat(Point2f(0.0f, kHeight))),
      1e-5f);
  EXPECT_VECTOR_NEAR(
      ion::math::Normalized(kRightBottomNear - eye_position_world),
      ion::math::Normalized(camera.RayDirectionFloat(Point2f(kWidth, 0.0f))),
      1e-5f);
  EXPECT_VECTOR_NEAR(
      ion::math::Normalized(kRightTopNear - eye_position_world),
      ion::math::Normalized(camera.RayDirectionFloat(Point2f(kWidth, kHeight))),
      1e-5f);

  // Test integer version.
  EXPECT_EQ(camera.RayDirectionFloat(Point2f(0.5f, 0.5f)),
            camera.RayDirection(Point2i(0, 0)));
  EXPECT_EQ(camera.RayDirectionFloat(Point2f(kWidth - 0.5f, kHeight - 0.5f)),
            camera.RayDirection(Point2i(kWidth - 1, kHeight - 1)));
}

TEST(ProjectiveCamera, RayEndWindowZ) {
  ProjectiveCamera camera =
      CreateTestCamera(ProjectiveCamera::DepthType::kWindowZ);

  // Verify that the 8 corner points of the frustum are handled correctly.
  EXPECT_VECTOR_NEAR(kLeftBottomNear,
                     camera.RayEndFloat(Point2f(0.0f, 0.0f), 0.0f),
                     2e-5f * kNear);
  EXPECT_VECTOR_NEAR(kLeftBottomFar,
                     camera.RayEndFloat(Point2f(0.0f, 0.0f), 1.0f),
                     2e-5f * kFar);
  EXPECT_VECTOR_NEAR(kLeftTopNear,
                     camera.RayEndFloat(Point2f(0.0f, kHeight), 0.0f),
                     2e-5f * kNear);
  EXPECT_VECTOR_NEAR(kLeftTopFar,
                     camera.RayEndFloat(Point2f(0.0f, kHeight), 1.0f),
                     2e-5f * kFar);
  EXPECT_VECTOR_NEAR(kRightBottomNear,
                     camera.RayEndFloat(Point2f(kWidth, 0.0f), 0.0f),
                     2e-5f * kNear);
  EXPECT_VECTOR_NEAR(kRightBottomFar,
                     camera.RayEndFloat(Point2f(kWidth, 0.0f), 1.0f),
                     2e-5f * kFar);
  EXPECT_VECTOR_NEAR(kRightTopNear,
                     camera.RayEndFloat(Point2f(kWidth, kHeight), 0.0f),
                     2e-5f * kNear);
  EXPECT_VECTOR_NEAR(kRightTopFar,
                     camera.RayEndFloat(Point2f(kWidth, kHeight), 1.0f),
                     2e-5f * kFar);

  // Test integer version.
  EXPECT_EQ(camera.RayEndFloat(Point2f(0.5f, 0.5f), 0.0f),
            camera.RayEnd(Point2i(0, 0), 0.0f));
  EXPECT_EQ(camera.RayEndFloat(Point2f(kWidth - 0.5f, kHeight - 0.5f), 1.0f),
            camera.RayEnd(Point2i(kWidth - 1, kHeight - 1), 1.0f));
}

TEST(ProjectiveCamera, RayEndEyeZ) {
  ProjectiveCamera camera =
      CreateTestCamera(ProjectiveCamera::DepthType::kEyeZ);

  // Verify that the 8 corner points of the frustum are handled correctly.
  EXPECT_VECTOR_NEAR(kLeftBottomNear,
                     camera.RayEndFloat(Point2f(0.0f, 0.0f), kNear), 2e-5f);
  EXPECT_VECTOR_NEAR(kLeftBottomFar,
                     camera.RayEndFloat(Point2f(0.0f, 0.0f), kFar), 2e-5f);
  EXPECT_VECTOR_NEAR(kLeftTopNear,
                     camera.RayEndFloat(Point2f(0.0f, kHeight), kNear), 2e-5f);
  EXPECT_VECTOR_NEAR(kLeftTopFar,
                     camera.RayEndFloat(Point2f(0.0f, kHeight), kFar), 2e-5f);
  EXPECT_VECTOR_NEAR(kRightBottomNear,
                     camera.RayEndFloat(Point2f(kWidth, 0.0f), kNear), 2e-5f);
  EXPECT_VECTOR_NEAR(kRightBottomFar,
                     camera.RayEndFloat(Point2f(kWidth, 0.0f), kFar), 2e-5f);
  EXPECT_VECTOR_NEAR(kRightTopNear,
                     camera.RayEndFloat(Point2f(kWidth, kHeight), kNear),
                     2e-5f);
  EXPECT_VECTOR_NEAR(kRightTopFar,
                     camera.RayEndFloat(Point2f(kWidth, kHeight), kFar), 2e-5f);

  // Test integer version.
  EXPECT_EQ(camera.RayEndFloat(Point2f(0.5f, 0.5f), kNear),
            camera.RayEnd(Point2i(0, 0), kNear));
  EXPECT_EQ(camera.RayEndFloat(Point2f(kWidth - 0.5f, kHeight - 0.5f), kFar),
            camera.RayEnd(Point2i(kWidth - 1, kHeight - 1), kFar));
}

TEST(ProjectiveCamera, RayEndRayDepth) {
  ProjectiveCamera camera =
      CreateTestCamera(ProjectiveCamera::DepthType::kRayDepth);

  // Verify that the 8 corner points of the frustum are handled correctly.
  const float kLeftBottomNearDepth =
      ion::math::Length(kLeftBottomNear - (Point3f::Zero() - kTranslation));
  const float kLeftBottomFarDepth =
      ion::math::Length(kLeftBottomFar - (Point3f::Zero() - kTranslation));
  const float kLeftTopNearDepth =
      ion::math::Length(kLeftTopNear - (Point3f::Zero() - kTranslation));
  const float kLeftTopFarDepth =
      ion::math::Length(kLeftTopFar - (Point3f::Zero() - kTranslation));
  const float kRightBottomNearDepth =
      ion::math::Length(kRightBottomNear - (Point3f::Zero() - kTranslation));
  const float kRightBottomFarDepth =
      ion::math::Length(kRightBottomFar - (Point3f::Zero() - kTranslation));
  const float kRightTopNearDepth =
      ion::math::Length(kRightTopNear - (Point3f::Zero() - kTranslation));
  const float kRightTopFarDepth =
      ion::math::Length(kRightTopFar - (Point3f::Zero() - kTranslation));

  EXPECT_VECTOR_NEAR(
      kLeftBottomNear,
      camera.RayEndFloat(Point2f(0.0f, 0.0f), kLeftBottomNearDepth), 3e-5f);
  EXPECT_VECTOR_NEAR(
      kLeftBottomFar,
      camera.RayEndFloat(Point2f(0.0f, 0.0f), kLeftBottomFarDepth), 3e-5f);
  EXPECT_VECTOR_NEAR(
      kLeftTopNear,
      camera.RayEndFloat(Point2f(0.0f, kHeight), kLeftTopNearDepth), 3e-5f);
  EXPECT_VECTOR_NEAR(
      kLeftTopFar, camera.RayEndFloat(Point2f(0.0f, kHeight), kLeftTopFarDepth),
      3e-5f);
  EXPECT_VECTOR_NEAR(
      kRightBottomNear,
      camera.RayEndFloat(Point2f(kWidth, 0.0f), kRightBottomNearDepth), 3e-5f);
  EXPECT_VECTOR_NEAR(
      kRightBottomFar,
      camera.RayEndFloat(Point2f(kWidth, 0.0f), kRightBottomFarDepth), 3e-5f);
  EXPECT_VECTOR_NEAR(
      kRightTopNear,
      camera.RayEndFloat(Point2f(kWidth, kHeight), kRightTopNearDepth), 3e-5f);
  EXPECT_VECTOR_NEAR(
      kRightTopFar,
      camera.RayEndFloat(Point2f(kWidth, kHeight), kRightTopFarDepth), 3e-5f);

  // Test integer version.
  EXPECT_EQ(camera.RayEndFloat(Point2f(0.5f, 0.5f), kLeftBottomNearDepth),
            camera.RayEnd(Point2i(0, 0), kLeftBottomNearDepth));
  EXPECT_EQ(camera.RayEndFloat(Point2f(kWidth - 0.5f, kHeight - 0.5f),
                               kRightTopFarDepth),
            camera.RayEnd(Point2i(kWidth - 1, kHeight - 1), kRightTopFarDepth));
}

TEST(ProjectiveCamera, NearFarOrthographic) {
  // Construct an orthographic camera.
  const ProjectiveCamera camera(kImageSize,
                                ion::math::OrthographicMatrixFromFrustum(
                                    kLeft, kRight, kBottom, kTop, kNear, kFar),
                                ion::math::Matrix4f::Identity());

  // Check that near and far clipping planes are computed correctly from the
  // orthographics clip-from-eye matrix.
  float near_clip, far_clip;
  camera.ComputeNearFar(&near_clip, &far_clip);
  EXPECT_NEAR(kNear, near_clip, 1e-3f);
  EXPECT_NEAR(kFar, far_clip, 1e-3f);
}

TEST(ProjectiveCamera, NearFarPerspective) {
  // Construct a perspective camera.
  const ProjectiveCamera camera(kImageSize,
                                ion::math::PerspectiveMatrixFromFrustum(
                                    kLeft, kRight, kBottom, kTop, kNear, kFar),
                                ion::math::Matrix4f::Identity());

  // Check that near and far clipping planes are computed correctly from the
  // perspective clip-from-eye matrix.
  float near_clip, far_clip;
  camera.ComputeNearFar(&near_clip, &far_clip);
  EXPECT_NEAR(kNear, near_clip, 1e-3f);
  EXPECT_NEAR(kFar, far_clip, 1e-3f);
}

TEST(ProjectiveCamera, ClipPlanesOrthographic) {
  // Construct an orthographic camera.
  const ProjectiveCamera camera(kImageSize,
                                ion::math::OrthographicMatrixFromFrustum(
                                    kLeft, kRight, kBottom, kTop, kNear, kFar),
                                ion::math::Matrix4f::Identity());

  // Check that the clipping planes are computed correctly from the
  // orthographics clip-from-eye matrix.
  float left, right, bottom, top, near_clip, far_clip;
  camera.ComputeClipPlanes(&left, &right, &bottom, &top, &near_clip, &far_clip);
  EXPECT_NEAR(kLeft, left, 1e-3f);
  EXPECT_NEAR(kRight, right, 1e-3f);
  EXPECT_NEAR(kBottom, bottom, 1e-3f);
  EXPECT_NEAR(kTop, top, 1e-3f);
  EXPECT_NEAR(kNear, near_clip, 1e-3f);
  EXPECT_NEAR(kFar, far_clip, 1e-3f);
}

TEST(ProjectiveCamera, ClipPlanesPerspective) {
  // Construct a perspective camera.
  const ProjectiveCamera camera(kImageSize,
                                ion::math::PerspectiveMatrixFromFrustum(
                                    kLeft, kRight, kBottom, kTop, kNear, kFar),
                                ion::math::Matrix4f::Identity());

  // Check that the clipping planes are computed correctly from the
  // perspective clip-from-eye matrix.
  float left, right, bottom, top, near_clip, far_clip;
  camera.ComputeClipPlanes(&left, &right, &bottom, &top, &near_clip, &far_clip);
  EXPECT_NEAR(kLeft, left, 1e-3f);
  EXPECT_NEAR(kRight, right, 1e-3f);
  EXPECT_NEAR(kBottom, bottom, 1e-3f);
  EXPECT_NEAR(kTop, top, 1e-3f);
  EXPECT_NEAR(kNear, near_clip, 1e-3f);
  EXPECT_NEAR(kFar, far_clip, 1e-3f);
}

}  // namespace base
}  // namespace seurat
