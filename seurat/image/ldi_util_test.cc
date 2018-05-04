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

#include "seurat/image/ldi_util.h"

#include <vector>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/color.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace image {
namespace {

using base::Color3f;
using base::Color4f;
using ion::math::Vector2i;

TEST(LdiUtil, FlattenLdi) {
  const Color4f kOpaqueGreen(0.0f, 1.0f, 0.0f, 1.0f);
  const Color4f kTransparentBlack(0.0f, 0.0f, 0.0f, 0.0f);
  const Color4f kSemiTransparentRed(0.5f, 0.0f, 0.0f, 0.5f);
  const Color4f kSemiTransparentBlue(0.0f, 0.0f, 0.5f, 0.5f);
  const Color4f kOpaqueRed(1.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kOpaqueBlue(0.0f, 0.0f, 1.0f, 1.0f);
  const Vector2i kSize(2, 3);
  const std::vector<int> kSampleCounts{{1, 2, 2, 2, 0, 1}};
  const std::vector<Color4f> kColors{
      {kOpaqueGreen,                       // pixel (0, 0)
       kOpaqueGreen, kSemiTransparentRed,  // pixel (1, 0)
       kSemiTransparentRed, kOpaqueBlue,   // pixel (0, 1)
       kSemiTransparentBlue, kOpaqueRed,   // pixel (1, 1)
                                           // pixel (0, 2) is empty
       kSemiTransparentRed}};              // pixel (1, 2)
  // Use arbitrary depth values. They should have no influence on the
  // result
  const std::vector<float> kDepths{
      {0.5f, 0.2f, 0.4f, 0.3f, 0.5f, 0.1f, 0.8f, 0.7f}};
  Ldi4f ldi(kSize, kSampleCounts, kColors, kDepths);
  Image3f image = FlattenLdi(ldi);
  EXPECT_EQ(kSize, image.GetSize());
  EXPECT_EQ(Color3f(0.0f, 1.0f, 0.0f), image.At(0, 0));
  EXPECT_EQ(Color3f(0.0f, 1.0f, 0.0f), image.At(1, 0));
  EXPECT_EQ(Color3f(0.5f, 0.0f, 0.5f), image.At(0, 1));
  EXPECT_EQ(Color3f(0.5f, 0.0f, 0.5f), image.At(1, 1));
  EXPECT_EQ(Color3f(0.0f, 0.0f, 0.0f), image.At(0, 2));
  EXPECT_EQ(Color3f(0.5f, 0.0f, 0.0f), image.At(1, 2));
}

TEST(LdiUtil, FlipLdiVertically) {
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kBlue(0.0f, 0.0f, 1.0f, 1.0f);
  const Vector2i kSize(2, 3);
  const std::vector<int> kSampleCounts{{1, 2, 2, 2, 0, 1}};
  const std::vector<Color4f> kColors{{kRed,          // pixel (0, 0)
                                      kRed, kBlue,   // pixel (1, 0)
                                      kRed, kRed,    // pixel (0, 1)
                                      kBlue, kBlue,  // pixel (1, 1)
                                                     // pixel (0, 2) is empty
                                      kBlue}};       // pixel (1, 2)
  // Use arbitrary depth values. They should have no influence on the
  // result
  const std::vector<float> kDepths{
      {0.5f, 0.2f, 0.4f, 0.3f, 0.5f, 0.1f, 0.8f, 0.7f}};
  Ldi4f ldi(kSize, kSampleCounts, kColors, kDepths);
  Ldi4f flipped = FlipLdiVertically(ldi);
  EXPECT_EQ(kSize, flipped.GetSize());
  EXPECT_EQ(0, flipped.GetSampleCount({0, 0}));

  EXPECT_EQ(1, flipped.GetSampleCount({1, 0}));
  EXPECT_EQ(kBlue, flipped.GetColors({1, 0})[0]);
  EXPECT_EQ(0.7f, flipped.GetDepths({1, 0})[0]);

  EXPECT_EQ(2, flipped.GetSampleCount({0, 1}));
  EXPECT_EQ(kRed, flipped.GetColors({0, 1})[0]);
  EXPECT_EQ(kRed, flipped.GetColors({0, 1})[1]);
  EXPECT_EQ(0.3f, flipped.GetDepths({0, 1})[0]);
  EXPECT_EQ(0.5f, flipped.GetDepths({0, 1})[1]);

  EXPECT_EQ(2, flipped.GetSampleCount({1, 1}));
  EXPECT_EQ(kBlue, flipped.GetColors({1, 1})[0]);
  EXPECT_EQ(kBlue, flipped.GetColors({1, 1})[1]);
  EXPECT_EQ(0.1f, flipped.GetDepths({1, 1})[0]);
  EXPECT_EQ(0.8f, flipped.GetDepths({1, 1})[1]);

  EXPECT_EQ(1, flipped.GetSampleCount({0, 2}));
  EXPECT_EQ(kRed, flipped.GetColors({0, 2})[0]);
  EXPECT_EQ(0.5f, flipped.GetDepths({0, 2})[0]);

  EXPECT_EQ(2, flipped.GetSampleCount({1, 2}));
  EXPECT_EQ(kRed, flipped.GetColors({1, 2})[0]);
  EXPECT_EQ(kBlue, flipped.GetColors({1, 2})[1]);
  EXPECT_EQ(0.2f, flipped.GetDepths({1, 2})[0]);
  EXPECT_EQ(0.4f, flipped.GetDepths({1, 2})[1]);
}

}  // namespace
}  // namespace image
}  // namespace seurat
