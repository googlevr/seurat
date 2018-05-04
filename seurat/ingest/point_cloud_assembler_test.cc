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

#include "seurat/ingest/point_cloud_assembler.h"

#include <memory>

#include "ion/math/matrix.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/base/color.h"
#include "seurat/base/projective_camera.h"
#include "seurat/geometry/binning_point_cloud_builder.h"
#include "seurat/geometry/point_cloud_builder.h"
#include "seurat/image/image.h"
#include "seurat/image/ldi.h"
#include "seurat/ingest/view_group_loader.h"
#include "seurat/ingest/view_group_loader_test_utils.h"

namespace seurat {
namespace ingest {
namespace {

using base::Color3f;
using base::Color4f;
using geometry::BinningPointCloudBuilder;
using geometry::PointCloudBuilder;
using image::Image1f;
using image::Image3f;
using ion::math::Matrix4f;
using ion::math::Point3f;
using ion::math::Vector2i;
using ion::math::Vector3f;
using ion::math::Vector3i;

TEST(PointCloudAssemblerTest, SmokeTest) {
  // This test runs the PointCloudAssembler with a FakeViewGroupLoader and
  // checks if the result has the expected number of points and cameras.

  const int kNumViewGroups = 1;
  const Vector2i kImageSize(16, 16);
  const std::array<Color3f, 6> kFaceColors = {{
      Color3f(0.0f, 0.0f, -1.0f),
      Color3f(0.0f, 0.0f, 1.0f),
      Color3f(-1.0f, 0.0f, 0.0f),
      Color3f(1.0f, 0.0f, 0.0f),
      Color3f(0.0f, -1.0f, 0.0f),
      Color3f(0.0f, 1.0f, 0.0f),
  }};
  std::array<float, 6> kFaceDepths = {{0.1f, 0.5f, 0.3f, 0.6f, 0.4f, 0.2f}};

  const int kThreadCount = 1;
  const Vector3i kResolution(16, 16, 256);
  const float kNear = 0.01f;
  std::unique_ptr<ViewGroupLoader> view_group_loader(new FakeViewGroupLoader(
      kNumViewGroups, kImageSize, kFaceColors, kFaceDepths));
  std::unique_ptr<PointCloudBuilder> point_cloud_builder(
      new BinningPointCloudBuilder(kThreadCount, kResolution, kNear));
  PointCloudAssembler point_cloud_assembler(kThreadCount,
                                            std::move(view_group_loader),
                                            std::move(point_cloud_builder));
  std::vector<Point3f> positions;
  std::vector<float> weights;
  EXPECT_TRUE(point_cloud_assembler.Build(&positions, &weights).ok());
  const int64 kExpectedNumberOfPoints =
      kImageSize[0] * kImageSize[1] * kNumViewGroups * 6;
  EXPECT_EQ(kExpectedNumberOfPoints, positions.size());
  for (float weight : weights) {
    EXPECT_EQ(1.0f, weight);
  }
}

TEST(PointCloudAssemblerTest, PositionsFromView) {
  const Vector2i kSize(2, 3);
  const float kLeft = -0.8f;
  const float kRight = 1.0f;
  const float kBottom = -0.7f;
  const float kTop = 0.9f;
  const float kNear = 1.0f;
  const float kFar = 100.0f;
  const Vector3f kTranslation(3.0f, -1.0f, 2.0f);
  const Matrix4f kClipFromEyeMatrix(ion::math::PerspectiveMatrixFromFrustum(
      kLeft, kRight, kBottom, kTop, kNear, kFar));
  const Matrix4f kEyeFromWorldMatrix(
      ion::math::TranslationMatrix(kTranslation));
  const base::ProjectiveCamera kCamera(kSize, kClipFromEyeMatrix,
                                       kEyeFromWorldMatrix);

  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  const std::vector<int> kSampleCounts{{1, 2, 2, 2, 0, 1}};
  const std::vector<Color4f> kColors(8, kRed);
  const std::vector<float> kDepths{{0.5f,        // pixel (0, 0)
                                    0.2f, 0.4f,  // pixel (1, 0)
                                    0.3f, 0.5f,  // pixel (0, 1)
                                    0.1f, 0.8f,  // pixel (1, 1)
                                                 // pixel (0, 2) is empty
                                    0.7f}};      // pixel (1, 2)
  const image::Ldi4f kLdi(kSize, kSampleCounts, kColors, kDepths);
  const int kThreadCount = 1;
  const std::vector<Point3f> actual_positions =
      PointCloudAssembler::PositionsFromView(kThreadCount, kCamera, kLdi);

  std::vector<Point3f> expected_positions{
      {kCamera.RayEnd({0, 0}, kDepths[0]), kCamera.RayEnd({1, 0}, kDepths[1]),
       kCamera.RayEnd({1, 0}, kDepths[2]), kCamera.RayEnd({0, 1}, kDepths[3]),
       kCamera.RayEnd({0, 1}, kDepths[4]), kCamera.RayEnd({1, 1}, kDepths[5]),
       kCamera.RayEnd({1, 1}, kDepths[6]), kCamera.RayEnd({1, 2}, kDepths[7])}};

  EXPECT_THAT(actual_positions,
              ::testing::UnorderedElementsAreArray(expected_positions));
}

}  // namespace
}  // namespace ingest
}  // namespace seurat
