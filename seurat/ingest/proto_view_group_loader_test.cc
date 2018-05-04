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

#include "seurat/ingest/proto_view_group_loader.h"

#include "ion/math/matrixutils.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/api/api.pb.h"
#include "seurat/api/math.pb.h"
#include "seurat/base/color.h"
#include "seurat/base/file_system.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/base/math_util.h"
#include "seurat/base/projective_camera.h"
#include "seurat/base/projective_camera_util.h"
#include "seurat/image/image_test_utils.h"
#include "seurat/image/ldi.h"
#include "seurat/image/ldi_test_utils.h"
#include "seurat/ingest/ldi_loader_test_utils.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace ingest {
namespace {

using base::Color4f;
using base::ProjectiveCamera;
using ion::math::Matrix4f;
using ion::math::Point3f;
using ion::math::Point4f;
using ion::math::Vector2i;
using ion::math::Vector3f;

void TestMatrices(const Matrix4f& actual_matrix,
                  const Matrix4f& expected_matrix, float matrix_epsilon) {
  for (int row = 0; row < 4; ++row) {
    EXPECT_VECTOR_NEAR(ion::math::Row(actual_matrix, row),
                       ion::math::Row(expected_matrix, row), matrix_epsilon);
  }
}

void TestProjectiveCameras(const ProjectiveCamera& actual_camera,
                           const ProjectiveCamera& expected_camera,
                           float matrix_epsilon) {
  EXPECT_EQ(actual_camera.GetImageSize(), expected_camera.GetImageSize());
  TestMatrices(actual_camera.GetEyeFromWorld(),
               expected_camera.GetEyeFromWorld(), matrix_epsilon);
  TestMatrices(actual_camera.GetClipFromEye(), expected_camera.GetClipFromEye(),
               matrix_epsilon);
}

TEST(ProtoViewGroupLoader, LoadEmptyCapture) {
  api::proto::Capture empty_capture;
  std::shared_ptr<base::FileSystem> file_system =
      std::make_shared<base::FileSystem>("/dev/null");
  std::unique_ptr<FakeLdiLoader> ldi_loader(new FakeLdiLoader(image::Ldi4f()));
  const int kSingleLoaderThread = 1;
  ProtoViewGroupLoader view_group_loader(empty_capture, kSingleLoaderThread,
                                         std::move(ldi_loader),
                                         std::move(file_system));
  EXPECT_EQ(0, view_group_loader.GetNumViewGroups());
  EXPECT_FALSE(view_group_loader.LoadViewGroup(0, nullptr, nullptr).ok());
}

TEST(ProtoViewGroupLoader, LoadWithFakeLdiLoader) {
  const int kNumViewGroups = 3;
  const int kNumViews = 5;
  const Vector2i kImageSize(2, 3);
  // Position the camera backward from the origin on Z, raised a bit.
  const Point3f kEyePositionWorld(0.0f, 0.1f, 2.0f);
  // Look toward -Z.
  const Point3f kEyeTargetWorld(0.0f, 0.0f, -2.0f);
  // Y-up, following Seurat convention.
  const Vector3f kYUpWorld(0.0f, 1.0f, 0.0f);
  const float kNearZ = 1.0f;
  const Matrix4f eye_from_world = ion::math::LookAtMatrixFromCenter(
      kEyePositionWorld, kEyeTargetWorld, kYUpWorld);
  const Matrix4f world_from_eye = ion::math::Inverse(eye_from_world);
  const ProjectiveCamera kCamera(
      kImageSize,
      base::PerspectiveMatrixFromInfiniteFrustum(-kNearZ, kNearZ, -kNearZ,
                                                 kNearZ, kNearZ, 0.0f),
      eye_from_world, world_from_eye, ProjectiveCamera::DepthType::kWindowZ);
  api::proto::Capture capture;
  for (int view_group_index = 0; view_group_index < kNumViewGroups;
       ++view_group_index) {
    api::proto::ViewGroup* view_group = capture.add_view_groups();
    for (int view_index = 0; view_index < kNumViews; ++view_index) {
      api::proto::View* view = view_group->add_views();
      base::ProjectiveCameraToProto(
          kCamera, view->mutable_camera()->mutable_projective());
    }
  }

  const std::vector<int> kSampleCounts(6, 1);
  const std::vector<Color4f> kColors(6, Color4f(1.0f, 0.0f, 0.0f, 1.0f));
  const std::vector<float> kDepths(6, 1.0f);
  const image::Ldi4f kLdi(kImageSize, kSampleCounts, kColors, kDepths);

  std::shared_ptr<base::FileSystem> file_system =
      std::make_shared<base::FileSystem>("/dev/null");
  std::unique_ptr<FakeLdiLoader> ldi_loader(new FakeLdiLoader(kLdi));
  const int kTwoLoaderThreads = 2;
  ProtoViewGroupLoader view_group_loader(capture, kTwoLoaderThreads,
                                         std::move(ldi_loader),
                                         std::move(file_system));

  EXPECT_EQ(kNumViewGroups, view_group_loader.GetNumViewGroups());
  for (int view_group_index = 0; view_group_index < kNumViewGroups;
       ++view_group_index) {
    std::vector<std::shared_ptr<base::Camera>> cameras;
    std::vector<image::Ldi4f> ldis;
    base::Status status =
        view_group_loader.LoadViewGroup(view_group_index, &cameras, &ldis);
    EXPECT_TRUE(status.ok());
    EXPECT_EQ(kNumViews, cameras.size());
    EXPECT_EQ(kNumViews, ldis.size());
    for (int view_index = 0; view_index < kNumViews; ++view_index) {
      std::shared_ptr<ProjectiveCamera> actual_camera =
          std::dynamic_pointer_cast<ProjectiveCamera>(cameras[view_index]);
      ASSERT_TRUE(actual_camera != nullptr);
      TestProjectiveCameras(*actual_camera, kCamera, 1.0e-6f);
      EXPECT_TRUE(image::LdiEquals(kLdi, ldis[view_index]));
      // All the camera pose matrices for a view group should share the same
      // translation.
      Point3f translation = ion::math::ProjectPoint(
          actual_camera->GetWorldFromEye(), Point3f(0.0f, 0.0f, 0.0f));
      EXPECT_VECTOR_NEAR(translation, kEyePositionWorld, 1.0e-5f);
    }
  }
}

TEST(ProtoViewGroupLoader, LoadWithoutHeadboxCenter) {
  const float kNearZ = 0.5f;
  const Matrix4f kClipFromEye = base::PerspectiveMatrixFromInfiniteFrustum(
      -kNearZ, kNearZ, -kNearZ, kNearZ, kNearZ, 0.0f);
  const Point3f kEyePositionWorld(1.0f, 2.0f, 3.0f);
  const Matrix4f kWorldFromEye =
      ion::math::TranslationMatrix(kEyePositionWorld);
  const Matrix4f kEyeFromWorld = ion::math::Inverse(kWorldFromEye);
  const Vector2i kImageSize(1, 1);
  const ProjectiveCamera kCamera(kImageSize, kClipFromEye, kEyeFromWorld,
                                 kWorldFromEye,
                                 ProjectiveCamera::DepthType::kRayDepth);
  api::proto::Capture capture;
  api::proto::ViewGroup* view_group = capture.add_view_groups();
  api::proto::View* view = view_group->add_views();
  base::ProjectiveCameraToProto(kCamera,
                                view->mutable_camera()->mutable_projective());

  std::shared_ptr<base::FileSystem> file_system =
      std::make_shared<base::FileSystem>("/dev/null");
  std::unique_ptr<FakeLdiLoader> ldi_loader(new FakeLdiLoader(image::Ldi4f()));
  const int kSingleLoaderThread = 1;
  ProtoViewGroupLoader view_group_loader(capture, kSingleLoaderThread,
                                         std::move(ldi_loader),
                                         std::move(file_system));

  EXPECT_EQ(1, view_group_loader.GetNumViewGroups());
  std::vector<std::shared_ptr<base::Camera>> cameras;
  EXPECT_TRUE(view_group_loader.LoadViewGroup(0, &cameras, nullptr).ok());
  EXPECT_EQ(1, cameras.size());
  EXPECT_EQ(kImageSize, cameras[0]->GetImageSize());
  TestMatrices(cameras[0]->GetWorldFromEye(),
               ion::math::TranslationMatrix(kEyePositionWorld), 1e-5f);
  EXPECT_VECTOR_NEAR(Point3f(0.0f, 0.0f, -kNearZ) + kEyePositionWorld,
                     cameras[0]->RayOriginFloat({0.5f, 0.5f}), 1e-5f);
  EXPECT_VECTOR_NEAR(Point3f(0.0f, 0.0f, -2.0f) + kEyePositionWorld,
                     cameras[0]->RayEndFloat({0.5f, 0.5f}, 2.0f), 1e-5f);
  EXPECT_VECTOR_NEAR(Vector3f(0.0f, 0.0f, -kNearZ),
                     cameras[0]->RayDirectionFloat({0.5f, 0.5f}), 1e-5f);
}

TEST(ProtoViewGroupLoader, LoadWithHeadboxCenter) {
  const float kNearZ = 0.5f;
  const Matrix4f kClipFromEye = base::PerspectiveMatrixFromInfiniteFrustum(
      -kNearZ, kNearZ, -kNearZ, kNearZ, kNearZ, 0.0f);
  const Point3f kEyePositionWorld(1.0f, 2.0f, 3.0f);
  const Matrix4f kWorldFromEye =
      ion::math::TranslationMatrix(kEyePositionWorld);
  const Matrix4f kEyeFromWorld = ion::math::Inverse(kWorldFromEye);
  const Vector2i kImageSize(1, 1);
  const ProjectiveCamera kCamera(kImageSize, kClipFromEye, kEyeFromWorld,
                                 kWorldFromEye,
                                 ProjectiveCamera::DepthType::kRayDepth);

  api::proto::Capture capture;
  base::Point3fToProto(kEyePositionWorld, capture.mutable_headbox_center());
  api::proto::ViewGroup* view_group = capture.add_view_groups();
  api::proto::View* view = view_group->add_views();
  base::ProjectiveCameraToProto(kCamera,
                                view->mutable_camera()->mutable_projective());

  std::shared_ptr<base::FileSystem> file_system =
      std::make_shared<base::FileSystem>("/dev/null");
  std::unique_ptr<FakeLdiLoader> ldi_loader(new FakeLdiLoader(image::Ldi4f()));
  const int kSingleLoaderThread = 1;
  ProtoViewGroupLoader view_group_loader(capture, kSingleLoaderThread,
                                         std::move(ldi_loader),
                                         std::move(file_system));

  EXPECT_EQ(1, view_group_loader.GetNumViewGroups());
  std::vector<std::shared_ptr<base::Camera>> cameras;
  EXPECT_TRUE(view_group_loader.LoadViewGroup(0, &cameras, nullptr).ok());
  EXPECT_EQ(1, cameras.size());
  EXPECT_EQ(kImageSize, cameras[0]->GetImageSize());
  TestMatrices(cameras[0]->GetWorldFromEye(), Matrix4f::Identity(), 1e-5f);
  EXPECT_VECTOR_NEAR(Point3f(0.0f, 0.0f, -kNearZ),
                     cameras[0]->RayOriginFloat({0.5f, 0.5f}), 1e-5f);
  EXPECT_VECTOR_NEAR(Point3f(0.0f, 0.0f, -2.0f),
                     cameras[0]->RayEndFloat({0.5f, 0.5f}, 2.0f), 1e-5f);
  EXPECT_VECTOR_NEAR(Vector3f(0.0f, 0.0f, -kNearZ),
                     cameras[0]->RayDirectionFloat({0.5f, 0.5f}), 1e-5f);
}

}  // namespace
}  // namespace ingest
}  // namespace seurat
