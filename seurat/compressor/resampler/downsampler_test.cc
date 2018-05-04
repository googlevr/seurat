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

#include <vector>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/compressor/resampler/box_downsampler.h"
#include "seurat/compressor/resampler/gl_linear_upsampler.h"
#include "seurat/image/image.h"
#include "seurat/image/image_test_utils.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace compressor {

using base::Color4f;
using image::Image4f;
using image::Image4ui8;
using ion::math::Vector2i;

namespace {

// Downsample a checkerboard image using a box filter and compare the result
// against expected pixel values.
TEST(BoxDownsampler, DownsampleSmall) {
  constexpr int kSquareSize = 2;
  const Vector2i kSourceSize = Vector2i(4, 4);
  constexpr int kTargetWidth = 3;
  constexpr int kTargetHeight = 3;
  const Vector2i kTargetSize = Vector2i(kTargetWidth, kTargetHeight);
  constexpr float kEpsilon = 1.0e-2f;
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kGreen(0.0f, 1.0f, 0.0f, 1.0f);

  const std::vector<Color4f> colors = {{kGreen, kRed}};
  Image4f source_image =
      image::MakeCheckerboard(kSourceSize, kSquareSize, colors);
  BoxDownsampler downsampler;
  Image4f target_image = downsampler.Resample(source_image, kTargetSize);
  const Color4f expected_image[kTargetHeight][kTargetWidth] = {
      {{0.0f, 1.0f, 0.0f, 1.0f},
       {0.5f, 0.5f, 0.0f, 1.0f},
       {1.0f, 0.0f, 0.0f, 1.0f}},
      {{0.5f, 0.5f, 0.0f, 1.0f},
       {0.5f, 0.5f, 0.0f, 1.0f},
       {0.5f, 0.5f, 0.0f, 1.0f}},
      {{1.0f, 0.0f, 0.0f, 1.0f},
       {0.5f, 0.5f, 0.0f, 1.0f},
       {0.0f, 1.0f, 0.0f, 1.0f}}};
  for (int y = 0; y < target_image.Height(); ++y) {
    for (int x = 0; x < target_image.Width(); ++x) {
      EXPECT_VECTOR_NEAR(expected_image[y][x], target_image.At(x, y), kEpsilon);
    }
  }
}

// Downsample a checkerboard image using a box filter and verify pixel
// values near centers of checkerboard squares.
TEST(BoxDownsampler, Downsample) {
  constexpr int kSquareSize = 4;
  const Vector2i kSourceSize = Vector2i(24, 28);
  const Vector2i kTargetSize = Vector2i(20, 24);
  constexpr float kEpsilon = 1.0e-5f;
  const Color4f kWhite(1.0f, 1.0f, 1.0f, 1.0f);
  const Color4f kBlack(0.0f, 0.0f, 0.0f, 1.0f);
  const std::vector<Color4f> colors = {{kWhite, kBlack}};
  Image4f source_image =
      image::MakeCheckerboard(kSourceSize, kSquareSize, colors);
  BoxDownsampler downsampler;
  Image4f target_image = downsampler.Resample(source_image, kTargetSize);

  Vector2i step;
  for (int i = 0; i < 2; ++i) {
    step[i] = static_cast<float>(kSourceSize[i] - 1) / kSquareSize *
              (kTargetSize[i] - 1);
  }
  for (int y = step[1] / 2; y < target_image.Height(); y += step[1]) {
    for (int x = step[0] / 2; x < target_image.Width(); x += step[0]) {
      Color4f expected_color = colors[(x / step[0]) + (y / step[1]) % 2];
      EXPECT_VECTOR_NEAR(expected_color, target_image.At(x, y), kEpsilon);
    }
  }
}

// Upsample and then downsample a checkerboard image. Compare the result against
// the initial image.
TEST(BoxDownsampler, DownsampleUpsampledSmall) {
  constexpr int kSquareSize = 1;
  const Vector2i kSourceSize = Vector2i(3, 5);
  const Vector2i kUpsampledSize = Vector2i(6, 10);
  constexpr float kEpsilon = 0.4f;
  const Color4f kWhite(1.0f, 1.0f, 1.0f, 1.0f);
  const Color4f kBlack(0.0f, 0.0f, 0.0f, 1.0f);
  const std::vector<Color4f> colors = {{kWhite, kBlack}};

  Image4f source_image =
      image::MakeCheckerboard(kSourceSize, kSquareSize, colors);
  GlLinearUpsampler upsampler;
  Image4f upsampled_image = upsampler.Resample(source_image, kUpsampledSize);
  BoxDownsampler downsampler;
  Image4f target_image = downsampler.Resample(upsampled_image, kSourceSize);

  for (int y = kSquareSize / 2; y < kSourceSize[1]; y += kSquareSize) {
    for (int x = kSquareSize / 2; x < kSourceSize[0]; x += kSquareSize) {
      EXPECT_VECTOR_NEAR(source_image.At(x, y), target_image.At(x, y), kEpsilon)
          << " x " << x << " y " << y;
    }
  }
}

// Upsample and then downsample a checkerboard image. Compare the result against
// the initial image.
TEST(BoxDownsampler, DownsampleUpsampled) {
  constexpr int kSquareSize = 4;
  const Vector2i kSourceSize = Vector2i(24, 32);
  const Vector2i kUpsampledSize = Vector2i(36, 48);
  constexpr float kEpsilon = 2.0e-3f;
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kGreen(0.0f, 1.0f, 0.0f, 1.0f);
  const std::vector<Color4f> colors = {{kRed, kGreen}};

  Image4f source_image =
      image::MakeCheckerboard(kSourceSize, kSquareSize, colors);
  GlLinearUpsampler upsampler;
  Image4f upsampled_image = upsampler.Resample(source_image, kUpsampledSize);
  BoxDownsampler downsampler;
  Image4f target_image = downsampler.Resample(upsampled_image, kSourceSize);
  for (int y = kSquareSize / 2; y < kSourceSize[1]; y += kSquareSize) {
    for (int x = kSquareSize / 2; x < kSourceSize[0]; x += kSquareSize) {
      EXPECT_VECTOR_NEAR(source_image.At(x, y), target_image.At(x, y), kEpsilon)
          << " x " << x << " y " << y;
    }
  }
}

// Downsample and then upsample a checkerboard image. Compare the result against
// the initial image.
TEST(BoxDownsampler, UpsampleDownsampled) {
  constexpr int kSquareSize = 4;
  const Vector2i kSourceSize = Vector2i(36, 48);
  const Vector2i kDownsampledSize = Vector2i(24, 32);
  constexpr float kEpsilon = 2.0e-1f;
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kGreen(0.0f, 1.0f, 0.0f, 1.0f);
  const std::vector<Color4f> colors = {{kRed, kGreen}};

  Image4f source_image =
      image::MakeCheckerboard(kSourceSize, kSquareSize, colors);
  BoxDownsampler downsampler;
  Image4f downsampled_image =
      downsampler.Resample(source_image, kDownsampledSize);
  GlLinearUpsampler upsampler;
  Image4f target_image = upsampler.Resample(downsampled_image, kSourceSize);
  for (int y = kSquareSize / 2; y < kSourceSize[1]; y += kSquareSize) {
    for (int x = kSquareSize / 2; x < kSourceSize[0]; x += kSquareSize) {
      EXPECT_VECTOR_NEAR(source_image.At(x, y), target_image.At(x, y), kEpsilon)
          << " x " << x << " y " << y;
    }
  }
}

}  // namespace

}  // namespace compressor
}  // namespace seurat
