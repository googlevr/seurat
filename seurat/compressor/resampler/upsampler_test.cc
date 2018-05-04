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

#include "seurat/compressor/resampler/gl_linear_upsampler.h"

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/image/image.h"
#include "seurat/image/image_test_utils.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace compressor {

using base::Color4f;
using image::Image4f;
using image::Image4ui8;
using ion::math::Point2i;
using ion::math::Vector2f;
using ion::math::Vector2i;

namespace {

TEST(GlLinearUpsampler, UpsampleConstantWhite) {
  const Vector2i kSourceSize(2, 2);
  const Vector2i kTargetSize(8, 8);
  const Color4f kWhite(1.0f, 1.0f, 1.0f, 1.0f);
  constexpr float kEpsilon = 1.0e-5f;

  Image4f source_image(kSourceSize, kWhite);

  GlLinearUpsampler upsampler;
  Image4f target_image = upsampler.Resample(source_image, kTargetSize);

  for (int y = 0; y < target_image.Height(); ++y) {
    for (int x = 0; x < target_image.Width(); ++x) {
      EXPECT_VECTOR_NEAR(kWhite, target_image.At(x, y), kEpsilon);
    }
  }
}

// Upsample a checkerboard image using the OpenGL GL_LINEAR interpolation
// algorithm compare the result against expected pixel values.
TEST(GlLinearUpsampler, UpsampleSmall) {
  constexpr int kSquareSize = 1;
  const Vector2i kSourceSize = Vector2i(2, 2);
  constexpr int kUpsampledWidth = 3;
  constexpr int kUpsampledHeight = 3;
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kGreen(0.0f, 1.0f, 0.0f, 1.0f);
  const Color4f kBrown(0.5f, 0.5f, 0.0f, 1.0f);
  constexpr float kEpsilon = 0.01f;

  Image4f source_image =
      image::MakeCheckerboard(kSourceSize, kSquareSize, {{kRed, kGreen}});
  GlLinearUpsampler upsampler;
  Image4f target_image = upsampler.Resample(
      source_image, Vector2i(kUpsampledWidth, kUpsampledHeight));

  std::vector<Point2i> expected_green_texels = {{2, 0}, {0, 2}};
  for (const auto& texel : expected_green_texels) {
    EXPECT_VECTOR_NEAR(kGreen, target_image.At(texel), kEpsilon) << texel;
  }

  std::vector<Point2i> expected_red_texels = {{0, 0}, {2, 2}};
  for (const auto& texel : expected_red_texels) {
    EXPECT_VECTOR_NEAR(kRed, target_image.At(texel), kEpsilon) << texel;
  }

  std::vector<Point2i> expected_brown_texels = {
      {1, 0}, {0, 1}, {1, 1}, {2, 1}, {1, 2}};
  for (const auto& texel : expected_brown_texels) {
    EXPECT_VECTOR_NEAR(kBrown, target_image.At(texel), kEpsilon) << texel;
  }
}

// Upsample a checkerboard image using the OpenGL GL_LINEAR interpolation
// algorithm compare the result against expected pixel values.
TEST(GlLinearUpsampler, HalfTexelBorder) {
  constexpr int kSquareSize = 1;
  const Vector2i kSourceSize = Vector2i(3, 5);
  constexpr int kUpsampledWidth = 12;
  constexpr int kUpsampledHeight = 20;
  const Color4f kWhite(1.0f, 1.0f, 1.0f, 1.0f);
  const Color4f kBlack(0.0f, 0.0f, 0.0f, 1.0f);

  Image4f source_image =
      image::MakeCheckerboard(kSourceSize, kSquareSize, {{kWhite, kBlack}});
  GlLinearUpsampler upsampler;
  Image4f target_image = upsampler.Resample(
      source_image, Vector2i(kUpsampledWidth, kUpsampledHeight));
}

// Upsample a checkerboard image using the OpenGL GL_LINEAR interpolation
// algorithm and verify pixel values near centers of checkerboard squares.

TEST(GlLinearUpsampler, Upsample) {
  constexpr int kSquareSize = 4;
  const Vector2i kSourceSize(20, 24);
  const Vector2i kTargetSize(24, 30);
  constexpr float kEpsilon = 2.0e-3f;
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kGreen(0.0f, 1.0f, 0.0f, 1.0f);
  std::vector<Color4f> colors = {{kRed, kGreen}};
  Image4f source_image =
      image::MakeCheckerboard(kSourceSize, kSquareSize, colors);
  GlLinearUpsampler upsampler;
  Image4f target_image = upsampler.Resample(source_image, kTargetSize);
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

}  // namespace

}  // namespace compressor
}  // namespace seurat
