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

#include "seurat/base/ion_util_no_gl.h"

#include "ion/math/range.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "gtest/gtest.h"

using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Point3d;
using ion::math::Point3f;
using ion::math::Point3i;
using ion::math::Point4f;
using ion::math::Range1i;
using ion::math::Range2f;
using ion::math::Range2i;
using ion::math::Range3d;
using ion::math::Range3f;
using ion::math::Range3i;
using ion::math::Vector3f;

namespace seurat {
namespace base {
namespace {

template <typename T>
void FillImage(ion::gfx::ImagePtr image, T value) {
  T* data = image->GetData()->GetMutableData<T>();
  std::fill(data, data + (image->GetDataSize() / sizeof(T)), value);
}

TEST(IonUtilNoGlTest, Point4FromPoint3) {
  const Point3f p(1.0f, 2.0f, 3.0f);
  const Point4f p_hom = Point4FromPoint3(p);
  EXPECT_EQ(p_hom[0], 1.0f);
  EXPECT_EQ(p_hom[1], 2.0f);
  EXPECT_EQ(p_hom[2], 3.0f);
  EXPECT_EQ(p_hom[3], 1.0f);
}

TEST(IonUtilNoGlTest, Point3FromPoint4) {
  const Point4f p_hom(2.0f, 4.0f, 8.0f, 2.0f);
  const Point3f p = Point3FromPoint4(p_hom);
  EXPECT_EQ(p[0], 1.0f);
  EXPECT_EQ(p[1], 2.0f);
  EXPECT_EQ(p[2], 4.0f);
}

TEST(IonUtilNoGlTest, Point3FromPoint4FailsWithZeroW) {
  const Point4f p_hom(2.0f, 4.0f, 8.0f, 0.0f);
  EXPECT_DEBUG_DEATH(Point3FromPoint4(p_hom), ".*");
}

TEST(IonUtilNoGlTest, ProjectAABB) {
  const Range3f aabb({-2.0f, -4.0f, -6.0f}, {3.0f, 6.0f, 9.0f});
  {
    const Vector3f offset(0.5f, -1.0f, 10.0f);
    EXPECT_EQ(Range3f(aabb.GetMinPoint() + offset, aabb.GetMaxPoint() + offset),
              ProjectAABB(ion::math::TranslationMatrix(offset), aabb));
  }
  {
    const float kEpsilon = 1e5f;
    const Range3f expected_aabb({-6.0f, -2.0f, -6.0f}, {4.0f, 3.0f, 9.0f});
    const Range3f rotated_aabb =
        ProjectAABB(ion::math::RotationMatrixAxisAngleH(
                        Vector3f(0.0f, 0.0f, 1.0f),
                        ion::math::Anglef::FromRadians(3.14159265359f * 0.5f)),
                    aabb);
    EXPECT_NEAR(0.0f,
                ion::math::Distance(expected_aabb.GetMinPoint(),
                                    rotated_aabb.GetMinPoint()),
                kEpsilon);
    EXPECT_NEAR(0.0f,
                ion::math::Distance(expected_aabb.GetMaxPoint(),
                                    rotated_aabb.GetMaxPoint()),
                kEpsilon);
  }
}

TEST(IonUtilNoGlTest, DiscreteFromContinuousPixel) {
  EXPECT_EQ(Point2i(-1, 2), DiscreteFromContinuousPixel(Point2f(-0.1f, 2.3f)));
  EXPECT_EQ(Point2i(0, 0), DiscreteFromContinuousPixel(Point2f(0.1f, 0.9f)));
  EXPECT_EQ(Point2i(0, 1), DiscreteFromContinuousPixel(Point2f(0.1f, 1.0f)));
  EXPECT_EQ(Point2i(-1, 1), DiscreteFromContinuousPixel(Point2f(-1.0f, 1.1f)));
}

TEST(IonUtilNoGlTest, CompareImagesEqual) {
  ion::gfx::ImagePtr image_rgb8i_16x16 =
      CreateImage(ion::gfx::Image::Format::kRgb8i, {16, 16});
  ion::gfx::ImagePtr image_rgb8i_8x8 =
      CreateImage(ion::gfx::Image::Format::kRgb8i, {8, 8});
  ion::gfx::ImagePtr image_rgba8i_8x8 =
      CreateImage(ion::gfx::Image::Format::kRgba8i, {8, 8});
  FillImage(image_rgb8i_16x16, 0);
  FillImage(image_rgb8i_8x8, 0);
  FillImage(image_rgba8i_8x8, 0);

  EXPECT_FALSE(CompareImagesEqual(image_rgb8i_8x8, image_rgba8i_8x8));
  EXPECT_FALSE(CompareImagesEqual(image_rgb8i_8x8, image_rgb8i_16x16));
  EXPECT_FALSE(
      CompareImagesEqual(image_rgb8i_8x8, ion::gfx::ImagePtr(nullptr)));
  EXPECT_FALSE(
      CompareImagesEqual(ion::gfx::ImagePtr(nullptr), image_rgb8i_8x8));

  EXPECT_TRUE(CompareImagesEqual(ion::gfx::ImagePtr(nullptr),
                                 ion::gfx::ImagePtr(nullptr)));
  EXPECT_TRUE(CompareImagesEqual(image_rgb8i_8x8, image_rgb8i_8x8));
  EXPECT_TRUE(CompareImagesEqual(image_rgb8i_16x16, image_rgb8i_16x16));

  ion::gfx::ImagePtr image_rgb8i_16x16_white =
      CreateImage(ion::gfx::Image::Format::kRgb8i, {16, 16});
  FillImage(image_rgb8i_16x16_white, int8(255));

  EXPECT_FALSE(CompareImagesEqual(image_rgb8i_16x16, image_rgb8i_16x16_white));
}

TEST(IonUtilNoGlTest, EnclosingPixelRange1i) {
  Range1i range(-1, 2);
  Range1i expected(-1, 2);
  Range1i actual = EnclosingPixelRange(range);
  EXPECT_EQ(expected, actual);
}

TEST(IonUtilNoGlTest, EnclosingPixelRange2f) {
  Range2f range(Point2f(-1.2f, 2.6f), Point2f(3.7f, 5.2f));
  Range2i expected(Point2i(-2, 2), Point2i(3, 5));
  Range2i actual = EnclosingPixelRange(range);
  EXPECT_EQ(expected, actual);
}

TEST(IonUtilNoGlTest, EnclosingPixelRange3d) {
  Range3d range(Point3d(-1.2, 2.6, 0.2), Point3d(3.7, 5.2, 8.5));
  Range3i expected(Point3i(-2, 2, 0), Point3i(3, 5, 8));
  Range3i actual = EnclosingPixelRange(range);
  EXPECT_EQ(expected, actual);
}

}  // namespace
}  // namespace base
}  // namespace seurat
