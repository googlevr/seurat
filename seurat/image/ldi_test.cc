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

#include "seurat/image/ldi.h"

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/color.h"
#include "seurat/image/ldi_test_utils.h"

namespace seurat {
namespace image {
namespace {

using ion::math::Point2i;
using ion::math::Vector2i;
using seurat::base::Color3f;

TEST(Ldi, Empty) {
  Vector2i size(3, 2);
  std::vector<int> sample_counts{{0, 0, 0, 0, 0, 0}};
  std::vector<Color3f> colors;
  std::vector<float> depths;
  Ldi3f ldi(size, sample_counts, colors, depths);

  EXPECT_TRUE(LdiValid(ldi));

  EXPECT_EQ(size, ldi.GetSize());
  EXPECT_EQ(size[0], ldi.GetWidth());
  EXPECT_EQ(size[1], ldi.GetHeight());
  EXPECT_EQ(0, ldi.GetSampleCount());
  EXPECT_TRUE(ldi.GetColors().empty());
  EXPECT_TRUE(ldi.GetMutableColors().empty());
  EXPECT_TRUE(ldi.GetDepths().empty());
  for (int y = 0; y < size[1]; ++y) {
    for (int x = 0; x < size[0]; ++x) {
      EXPECT_EQ(0, ldi.GetSampleCount({x, y}));
      EXPECT_TRUE(ldi.GetColors({x, y}).empty());
      EXPECT_TRUE(ldi.GetMutableColors({x, y}).empty());
      EXPECT_TRUE(ldi.GetDepths({x, y}).empty());
    }
  }
}

TEST(Ldi, OneSamplePerPixel) {
  Vector2i size(3, 2);
  std::vector<Color3f> colors{
      {Color3f(0.1f, 0.1f, 0.1f), Color3f(0.2f, 0.2f, 0.2f),
       Color3f(0.3f, 0.3f, 0.3f), Color3f(0.4f, 0.4f, 0.4f),
       Color3f(0.5f, 0.5f, 0.5f), Color3f(0.6f, 0.6f, 0.6f)}};
  std::vector<float> depths{{0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f}};
  Ldi3f ldi(size, colors, depths);

  EXPECT_TRUE(LdiValid(ldi));

  EXPECT_EQ(size, ldi.GetSize());
  EXPECT_EQ(size[0], ldi.GetWidth());
  EXPECT_EQ(size[1], ldi.GetHeight());
  EXPECT_EQ(6, ldi.GetSampleCount());
  for (int y = 0; y < size[1]; ++y) {
    for (int x = 0; x < size[0]; ++x) {
      EXPECT_EQ(1, ldi.GetSampleCount({x, y}));
      EXPECT_EQ(colors[y * size[0] + x], ldi.GetColors({x, y})[0]);
      EXPECT_EQ(depths[y * size[0] + x], ldi.GetDepths({x, y})[0]);
    }
  }
}

TEST(Ldi, AccessColorsAndDepth) {
  // Create a test LDI with the following number of samples:
  // 2 0 2
  // 1 1 0
  //
  // First samples are red at depth=0.1, and second samples are blue at
  // depth=0.2.
  Vector2i size(3, 2);
  std::vector<int> sample_counts{{2, 0, 2, 1, 1, 0}};
  const Color3f kRed(1.0f, 0.0f, 0.0f);
  const Color3f kBlue(0.0f, 0.0f, 1.0f);
  std::vector<Color3f> colors{kRed, kBlue, kRed, kBlue, kRed, kRed};
  std::vector<float> depths{{0.1f, 0.2f, 0.1f, 0.2f, 0.1f, 0.1f}};
  Ldi3f ldi(size, sample_counts, colors, depths);

  EXPECT_TRUE(LdiValid(ldi));

  EXPECT_EQ(size, ldi.GetSize());
  EXPECT_EQ(size[0], ldi.GetWidth());
  EXPECT_EQ(size[1], ldi.GetHeight());
  EXPECT_EQ(6, ldi.GetSampleCount());
  EXPECT_EQ(absl::Span<const Color3f>(colors), ldi.GetColors());
  EXPECT_EQ(absl::Span<const float>(depths), ldi.GetDepths());
  for (int y = 0; y < size[1]; ++y) {
    for (int x = 0; x < size[0]; ++x) {
      const int index = y * size[0] + x;
      const int sample_count = sample_counts[index];
      EXPECT_EQ(sample_count, ldi.GetSampleCount({x, y}));
      EXPECT_EQ(sample_count, ldi.GetColors({x, y}).size());
      EXPECT_EQ(sample_count, ldi.GetDepths({x, y}).size());
      for (int s = 0; s < sample_count; ++s) {
        if (s == 0) {
          EXPECT_EQ(kRed, ldi.GetColors({x, y})[s]);
          EXPECT_EQ(0.1f, ldi.GetDepths({x, y})[s]);
        } else if (s == 1) {
          EXPECT_EQ(kBlue, ldi.GetColors({x, y})[s]);
          EXPECT_EQ(0.2f, ldi.GetDepths({x, y})[s]);
        } else {
          ASSERT_TRUE(false);
        }
      }
    }
  }
}

TEST(Ldi, MutateColors) {
  // Create a test LDI with the following number of samples:
  // 2 0 2
  // 1 1 0
  //
  // First samples are red at depth=0.1, and second samples are blue at
  // depth=0.2.
  Vector2i size(3, 2);
  std::vector<int> sample_counts{{2, 0, 2, 1, 1, 0}};
  const Color3f kRed(1.0f, 0.0f, 0.0f);
  const Color3f kBlue(0.0f, 0.0f, 1.0f);
  std::vector<Color3f> colors{kRed, kBlue, kRed, kBlue, kRed, kRed};
  std::vector<float> depths{{0.1f, 0.2f, 0.1f, 0.2f, 0.1f, 0.1f}};
  Ldi3f ldi(size, sample_counts, colors, depths);

  EXPECT_TRUE(LdiValid(ldi));

  absl::Span<Color3f> mutable_colors = ldi.GetMutableColors();
  for (auto& color : mutable_colors) {
    color *= 2.0f;
  }
  for (int i = 0; i < colors.size(); ++i) {
    EXPECT_EQ(colors[i] * 2.0f, ldi.GetColors()[i]);
  }
  absl::Span<Color3f> pixel_0_1_colors = ldi.GetMutableColors({0, 1});
  for (auto& color : pixel_0_1_colors) {
    color *= 2.0f;
  }
  for (int i = 0; i < colors.size(); ++i) {
    if (i == 4) {
      EXPECT_EQ(colors[i] * 4.0f, ldi.GetColors()[i]);
    } else {
      EXPECT_EQ(colors[i] * 2.0f, ldi.GetColors()[i]);
    }
  }
}

}  // namespace
}  // namespace image
}  // namespace seurat
