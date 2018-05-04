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

#include "seurat/ingest/clamping_view_group_loader.h"

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
using image::Ldi4f;
using ion::math::Point3f;
using ion::math::Vector2i;

void ExpectAllPointsWithinCube(const ViewGroupLoader& loader,
                               float cube_radius) {
  std::vector<std::shared_ptr<Camera>> cameras;
  std::vector<Ldi4f> ldis;
  for (int view_group = 0; view_group < loader.GetNumViewGroups();
       ++view_group) {
    EXPECT_TRUE(loader.LoadViewGroup(view_group, &cameras, &ldis).ok());

    for (int view_index = 0; view_index < cameras.size(); ++view_index) {
      const std::shared_ptr<Camera>& camera = cameras[view_index];
      const Ldi4f& ldi = ldis[view_index];
      for (int y = 0; y < ldi.GetHeight(); ++y) {
        for (int x = 0; x < ldi.GetWidth(); ++x) {
          for (float d : ldi.GetDepths({x, y})) {
            Point3f point = camera->RayEnd({x, y}, d);
            EXPECT_GT(cube_radius, std::fabs(point[0]));
            EXPECT_GT(cube_radius, std::fabs(point[1]));
            EXPECT_GT(cube_radius, std::fabs(point[2]));
          }
        }
      }
    }
  }
}

TEST(ClampingViewGroupLoaderTest, LoadsViews) {
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
    ClampingViewGroupLoader clamping_loader(
        4.0f, false,
        std::unique_ptr<ViewGroupLoader>(new FakeViewGroupLoader(
            kNumViewGroups, kImageSize, kFaceColors, kFaceDepths)));
    const float kEpsilon = 1.0e-4f;
    ExpectAllPointsWithinCube(clamping_loader, 4.0f + kEpsilon);
  }
  {
    ClampingViewGroupLoader clamping_loader_zero_is_inf(
        4.0f, true,
        std::unique_ptr<ViewGroupLoader>(new FakeViewGroupLoader(
            kNumViewGroups, kImageSize, kFaceColors, kFaceDepths)));
    const float kEpsilon = 1.0e-4f;
    ExpectAllPointsWithinCube(clamping_loader_zero_is_inf, 4.0f + kEpsilon);
  }
}

}  // namespace
}  // namespace ingest
}  // namespace seurat
