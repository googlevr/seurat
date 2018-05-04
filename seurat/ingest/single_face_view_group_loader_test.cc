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

#include "seurat/ingest/single_face_view_group_loader.h"

#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/base/camera.h"
#include "seurat/base/color.h"
#include "seurat/base/ion_util_no_gl.h"
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

TEST(SingleFaceViewGroupLoaderTest, LoadsViews) {
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

  const int kThreadCount = 3;

  SingleFaceViewGroupLoader single_face_loader(
      geometry::CubeFace::kLeft, kThreadCount,
      std::unique_ptr<ViewGroupLoader>(new FakeViewGroupLoader(
          kNumViewGroups, kImageSize, kFaceColors, kFaceDepths)));

  std::vector<std::shared_ptr<Camera>> cameras;
  std::vector<Ldi4f> ldis;
  EXPECT_TRUE(single_face_loader.LoadViewGroup(1, &cameras, &ldis).ok());
  EXPECT_EQ(ldis.size(), cameras.size());
  for (int i = 0; i < ldis.size(); ++i) {
    const auto& camera = cameras.at(i);
    const Ldi4f& ldi = ldis.at(i);

    for (int y = 0; y < ldi.GetHeight(); ++y) {
      for (int x = 0; x < ldi.GetWidth(); ++x) {
        for (int s = 0; s < ldi.GetSampleCount({x, y}); ++s) {
          Point3f p = camera->RayEnd({x, y}, ldi.GetDepths({x, y})[s]);
          // Ensure this is the left cube face, so x-coordinates are negative.
          EXPECT_GT(0.0f, p[0]);
          // The x-axis should be the major axis for all points.
          EXPECT_EQ(0, base::MajorAxisFromPosition(p));
        }
      }
    }
  }
}

}  // namespace
}  // namespace ingest
}  // namespace seurat
