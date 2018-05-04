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

#include "seurat/compressor/image_metrics/image_metrics.h"

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/image/image.h"
#include "seurat/image/image_test_utils.h"

namespace seurat {
namespace compressor {

using base::Color1f;
using base::Color4f;
using image::Image1f;
using image::Image3f;
using image::Image4f;
using image::Image4ui8;
using ion::math::Point2i;
using ion::math::Range2i;
using ion::math::Vector2i;

namespace {

TEST(ImageMetrics, EvalBitrate) {
  constexpr int kWidth = 12;
  constexpr int kHeight = 24;
  constexpr float kTolerance = 1.0e-5f;
  const Vector2i image_size(kWidth, kHeight);
  Image4f image(image_size, base::Color4f(0.1f, 0.3f, 0.9f, 1.0f));
  EXPECT_NEAR(9216.0f,
              EvalBitrate(image_size, ion::gfx::Image::Format::kRgba8888),
              kTolerance);
}

TEST(ImageMetrics, EvalRangeSSD) {
  const Vector2i kSize(5, 8);
  constexpr int kSquareSize = 2;
  constexpr float kTolerance = 1.0e-5f;
  const Color4f kWhite(1.0f, 1.0f, 1.0f, 1.0f);
  const Color4f kLightGray(0.9f, 0.9f, 0.9f, 1.0f);
  Image4f solid_white(kSize, kWhite);
  Image4f board =
      image::MakeCheckerboard(kSize, kSquareSize, {{kWhite, kLightGray}});
  Range2i range(Point2i(1, 2), Point2i(4, 6));
  EXPECT_NEAR(0.3f, EvalRangeSSD(solid_white, board, range), kTolerance);
}

TEST(ImageMetrics, EvalSSD) {
  const Vector2i kSize(5, 8);
  constexpr int kSquareSize = 2;
  constexpr float kTolerance = 1.0e-5f;
  const Color4f kWhite(1.0f, 1.0f, 1.0f, 1.0f);
  const Color4f kLightGray(0.9f, 0.9f, 0.9f, 1.0f);
  Image4f solid_white(kSize, kWhite);
  Image4f board =
      image::MakeCheckerboard(kSize, kSquareSize, {{kWhite, kLightGray}});
  EXPECT_NEAR(0.6f, EvalSSD(solid_white, board), kTolerance);
}

TEST(ImageMetrics, EvalPSNR) {
  const Vector2i kSize(5, 8);
  constexpr int kSquareSize = 2;
  constexpr float kTolerance = 1.0e-2f;
  const Color4f kWhite(1.0f, 1.0f, 1.0f, 1.0f);
  const Color4f kLightGray(0.9f, 0.9f, 0.9f, 1.0f);
  Image4f solid_white(kSize, kWhite);
  Image4f board =
      image::MakeCheckerboard(kSize, kSquareSize, {{kWhite, kLightGray}});
  EXPECT_NEAR(24.26f, EvalPSNR(solid_white, board), kTolerance);
}
}  // namespace

}  // namespace compressor
}  // namespace seurat
