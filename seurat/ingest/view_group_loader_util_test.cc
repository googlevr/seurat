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

#include "seurat/ingest/view_group_loader_util.h"

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/camera.h"
#include "seurat/base/color.h"
#include "seurat/ingest/view_group_loader_test_utils.h"

namespace seurat {
namespace ingest {
namespace {

using base::Camera;
using base::Color3f;
using image::Ldi4f;
using ion::math::Point3f;
using ion::math::Vector2i;

TEST(ViewGroupLoaderUtil, ComputeHeadbox) {
  const std::vector<ion::math::Point3f> kPositions{{Point3f(-0.1, -0.3f, -0.2f),
                                                    Point3f(0.2f, 0.3f, 0.5f),
                                                    Point3f(0.0f, 0.0f, 0.0f)}};
  const ion::math::Vector2i kImageSize(1, 1);
  FakeViewGroupLoader view_group_loader(kPositions, kImageSize);
  ion::math::Range3f headbox = ComputeHeadbox(view_group_loader);
  EXPECT_EQ(ion::math::Range3f({-0.1f, -0.3f, -0.2f}, {0.2f, 0.3f, 0.5f}),
            headbox);
}

TEST(ViewGroupLoaderUtil, ForEachViewGroupPrefetching) {
  const int kNumViewGroups = 20;
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

  FakeViewGroupLoader loader(kNumViewGroups, kImageSize, kFaceColors,
                             kFaceDepths);

  int counter = 0;
  base::Status status = ForEachViewGroupPrefetching(
      loader, [&](const std::vector<std::shared_ptr<Camera>>& cameras,
                  const std::vector<Ldi4f>& ldis) {
        counter++;
        return base::OkStatus();
      });
  EXPECT_TRUE(status.ok());
  EXPECT_EQ(counter, kNumViewGroups);
}

}  // namespace
}  // namespace ingest
}  // namespace seurat
