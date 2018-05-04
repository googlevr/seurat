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

#include "seurat/ingest/sphere_view_group_loader.h"

#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/base/camera.h"
#include "seurat/base/color.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/base/projective_camera.h"
#include "seurat/image/ldi.h"
#include "seurat/ingest/view_group_loader_test_utils.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace ingest {
namespace {

using base::Camera;
using base::Color3f;
using base::Color4f;
using image::Ldi4f;
using ion::math::Point3f;
using ion::math::Vector2i;

TEST(SphereViewGroupLoaderTest, NoSpheres) {
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

  SphereViewGroupLoader::SphereShader shader =
      [](const SphereViewGroupLoader::Sphere& sphere,
         const ion::math::Point3f& hit_point) -> Color3f {
    return Color3f::Zero();
  };

  // Create a SphereViewGroupLoader to test which does not actually render
  // any spheres.
  SphereViewGroupLoader sphere_loader(
      kThreadCount,
      std::unique_ptr<ViewGroupLoader>(new FakeViewGroupLoader(
          kNumViewGroups, kImageSize, kFaceColors, kFaceDepths)),
      std::vector<SphereViewGroupLoader::Sphere>{}, shader);

  // Since there are no spheres, the result of sphere_loader should be the
  // same as with expected_loader.
  //
  // Note that the depth values & camera-depth parameterization may (will!)
  // change, but the impact should be a noop.
  FakeViewGroupLoader expected_loader(kNumViewGroups, kImageSize, kFaceColors,
                                      kFaceDepths);

  std::vector<std::shared_ptr<Camera>> cameras;
  std::vector<Ldi4f> ldis;
  EXPECT_TRUE(sphere_loader.LoadViewGroup(1, &cameras, &ldis).ok());

  std::vector<std::shared_ptr<Camera>> expected_cameras;
  std::vector<Ldi4f> expected_ldis;
  ASSERT_TRUE(
      expected_loader.LoadViewGroup(1, &expected_cameras, &expected_ldis).ok());

  for (int i = 0; i < ldis.size(); ++i) {
    const auto& camera = cameras.at(i);
    const Ldi4f& ldi = ldis.at(i);
    const auto& expected_camera = expected_cameras.at(i);
    const Ldi4f& expected_ldi = expected_ldis.at(i);

    for (int y = 0; y < ldi.GetHeight(); ++y) {
      for (int x = 0; x < ldi.GetWidth(); ++x) {
        for (int s = 0; s < ldi.GetSampleCount({x, y}); ++s) {
          // Compare the ray endpoints.
          //
          // This allows the test to pass even if the ldi depths are
          // different, so long as the resulting camera compensates for the
          // change.
          Point3f point = camera->RayEnd({x, y}, ldi.GetDepths({x, y})[s]);
          Point3f expected_point = expected_camera->RayEnd(
              {x, y}, expected_ldi.GetDepths({x, y})[s]);
          EXPECT_VECTOR_NEAR(expected_point, point, 1.0e-5f);
          Color4f color = ldi.GetColors({x, y})[s];
          Color4f expected_color = expected_ldi.GetColors({x, y})[s];
          EXPECT_VECTOR_NEAR(expected_color, color, 1.0e-5f);
        }
      }
    }
  }
}

TEST(SphereViewGroupLoaderTest, RenderSeveralSpheres) {
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

  SphereViewGroupLoader::SphereShader shader =
      [](const SphereViewGroupLoader::Sphere& sphere,
         const ion::math::Point3f& hit_point) -> Color3f {
    return Color3f(1.0f, 0.0f, 0.0f);
  };

  std::vector<SphereViewGroupLoader::Sphere> spheres;
  spheres.push_back({{0.0f, 0.0f, -2.0f}, 2.0f});
  spheres.push_back({{0.0f, 0.0f, -3.0f}, 2.0f});
  spheres.push_back({{1.0f, 0.0f, 5.0f}, 5.0f});
  SphereViewGroupLoader sphere_loader(
      kThreadCount,
      std::unique_ptr<ViewGroupLoader>(new FakeViewGroupLoader(
          kNumViewGroups, kImageSize, kFaceColors, kFaceDepths)),
      std::move(spheres), shader);

  std::vector<std::shared_ptr<Camera>> cameras;
  std::vector<Ldi4f> ldis;
  // At least verify that it does not crash.
  //
  // Validating the rendered results look correct is out of scope for this test.
  EXPECT_TRUE(sphere_loader.LoadViewGroup(1, &cameras, &ldis).ok());
}

}  // namespace
}  // namespace ingest
}  // namespace seurat
