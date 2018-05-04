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

#include "seurat/ingest/subset_view_group_loader.h"

#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/base/camera.h"
#include "seurat/base/color.h"
#include "seurat/base/projective_camera.h"
#include "seurat/image/ldi.h"
#include "seurat/ingest/view_group_loader_test_utils.h"

namespace seurat {
namespace ingest {
namespace {

using base::Camera;
using base::Color3f;
using base::ProjectiveCamera;
using image::Ldi4f;
using ion::math::Vector2i;

TEST(SubsetViewGroupLoaderTest, LoadsViews) {
  const int kNumViewGroups = 3;
  const Vector2i kImageSize(3, 4);
  const std::array<Color3f, 6> kFaceColors = {{
      {1.0f, 0.0f, 0.0f},  //
      {0.0f, 1.0f, 0.0f},  //
      {0.0f, 0.0f, 1.0f},  //
      {1.0f, 1.0f, 0.0f},  //
      {0.0f, 1.0f, 1.0f},  //
      {0.0f, 0.0f, 0.0f}   //
  }};
  const std::array<float, 6> kFaceDepths = {
      {0.1f, 0.5f, 0.3f, 0.6f, 0.4f, 0.2f}};

  {
    SubsetViewGroupLoader single_loader(
        std::unique_ptr<ViewGroupLoader>(new FakeViewGroupLoader(
            kNumViewGroups, kImageSize, kFaceColors, kFaceDepths)),
        1);
    std::vector<std::shared_ptr<Camera>> cameras;
    std::vector<Ldi4f> ldis;
    EXPECT_TRUE(single_loader.LoadViewGroup(0, &cameras, &ldis).ok());
    EXPECT_EQ(6, cameras.size());
    EXPECT_EQ(6, ldis.size());
    EXPECT_FALSE(single_loader.LoadViewGroup(1, &cameras, &ldis).ok());
  }

  {
    std::unique_ptr<ViewGroupLoader> loader(new FakeViewGroupLoader(
        kNumViewGroups, kImageSize, kFaceColors, kFaceDepths));
    std::vector<std::shared_ptr<Camera>> expected_cameras[kNumViewGroups];
    for (int view_group_index = 0; view_group_index < kNumViewGroups;
         ++view_group_index) {
      base::Status status = loader->LoadViewGroup(
          view_group_index, &expected_cameras[view_group_index], nullptr);
      EXPECT_TRUE(status.ok());
    }
    SubsetViewGroupLoader all_loader(std::move(loader),
                                     std::numeric_limits<int>::max());
    EXPECT_EQ(3, all_loader.GetNumViewGroups());
    std::vector<std::shared_ptr<Camera>> actual_cameras;
    std::vector<Ldi4f> actual_ldis;
    for (int view_group_index = 0; view_group_index < kNumViewGroups;
         ++view_group_index) {
      base::Status status = all_loader.LoadViewGroup(
          view_group_index, &actual_cameras, &actual_ldis);
      EXPECT_TRUE(status.ok());
      EXPECT_EQ(6, actual_cameras.size());
      EXPECT_EQ(6, actual_ldis.size());
      for (int view_index = 0; view_index < 6; ++view_index) {
        std::shared_ptr<ProjectiveCamera> expected_camera =
            std::dynamic_pointer_cast<ProjectiveCamera>(
                expected_cameras[view_group_index][view_index]);
        std::shared_ptr<ProjectiveCamera> actual_camera =
            std::dynamic_pointer_cast<ProjectiveCamera>(
                actual_cameras[view_index]);
        ASSERT_TRUE(expected_camera != nullptr);
        ASSERT_TRUE(actual_camera != nullptr);
        EXPECT_EQ(*expected_camera, *actual_camera);
      }
    }
    EXPECT_FALSE(
        all_loader.LoadViewGroup(3, &actual_cameras, &actual_ldis).ok());
  }
}

}  // namespace
}  // namespace ingest
}  // namespace seurat
