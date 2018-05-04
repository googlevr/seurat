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

#include "seurat/compressor/rgba/rgba_encoding_selector.h"

#include <memory>
#include <vector>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/parallel.h"
#include "seurat/base/util.h"
#include "seurat/compressor/image_metrics/image_metrics.h"
#include "seurat/compressor/resampler/box_downsampler.h"
#include "seurat/compressor/resampler/gl_linear_upsampler.h"
#include "seurat/compressor/rgba/rgba_rate_resizer.h"
#include "seurat/image/atlaser.h"
#include "seurat/image/fixed_width_atlaser.h"
#include "seurat/image/image.h"
#include "seurat/image/image_test_utils.h"
#include "seurat/image/image_util.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace compressor {

using base::Color4f;
using image::FixedWidthAtlaser;
using image::Image4f;
using ion::math::Vector2i;

namespace {

// Generates a set of texture encodings.
std::vector<std::vector<EncodedTexture>> GenerateEncodings() {
  constexpr int kSquareSize = 2;
  const int kThreadCount = base::GetNumberOfHardwareThreads();
  const Vector2i kSourceSize = Vector2i(9, 12);
  const Vector2i kBlockSize(1, 1);
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kGreen(0.0f, 1.0f, 0.0f, 1.0f);
  const Color4f kBlue(0.0f, 0.0f, 1.0f, 1.0f);
  const std::vector<Color4f> colors = {{kRed, kGreen, kBlue}};
  Image4f texture = image::MakeCheckerboard(kSourceSize, kSquareSize, colors);
  RgbaRateResizer encoder(
      kThreadCount, std::unique_ptr<Resampler>(new GlLinearUpsampler()),
      std::unique_ptr<Resampler>(new BoxDownsampler()), kBlockSize);
  std::vector<std::vector<EncodedTexture>> resizings(1);
  encoder.Encode({texture}, absl::MakeSpan(resizings));
  return resizings;
}

// Sets a low target for the total bitrate. A low resolution version of the
// texture should be selected.
TEST(RgbaEncodingSelector, SmallAtlasSizeTarget) {
  const Vector2i kAtlasSizeTarget(2, 2);
  auto encodings = GenerateEncodings();
  auto atlaser = std::make_shared<FixedWidthAtlaser>(kAtlasSizeTarget);
  AtlasSizeTargetSelector encoding_selector(atlaser);
  std::vector<EncodedTexture> encoded_textures(encodings.size());
  encoding_selector.Select(encodings, absl::MakeSpan(encoded_textures));
  EXPECT_EQ(1, encoded_textures.size());
  EXPECT_GE(kAtlasSizeTarget[0], encoded_textures[0].image.Width());
  EXPECT_GE(kAtlasSizeTarget[1], encoded_textures[0].image.Height());
  EXPECT_NEAR(128.0f, encoded_textures[0].rate, 1.0e-3f);
  EXPECT_NEAR(72.2f, encoded_textures[0].distortion, 0.1f);
}

// Sets a constraint on the size of the selected texture.
TEST(RgbaEncodingSelector, MediumAtlasSizeTarget) {
  const Vector2i kAtlasSizeTarget(6, 9);
  auto encodings = GenerateEncodings();
  auto atlaser = std::make_shared<FixedWidthAtlaser>(kAtlasSizeTarget);
  AtlasSizeTargetSelector encoding_selector(atlaser);
  std::vector<EncodedTexture> encoded_textures(encodings.size());
  encoding_selector.Select(encodings, absl::MakeSpan(encoded_textures));
  EXPECT_EQ(1, encoded_textures.size());
  EXPECT_GE(kAtlasSizeTarget[0], encoded_textures[0].image.Width());
  EXPECT_GE(kAtlasSizeTarget[1], encoded_textures[0].image.Height());
  EXPECT_NEAR(960.0f, encoded_textures[0].rate, 1.0e-3f);
  EXPECT_NEAR(39.7f, encoded_textures[0].distortion, 0.1f);
}

// Sets a constraint on the size of the selected texture.
TEST(RgbaEncodingSelector, LargeAtlasSizeTarget) {
  const Vector2i kAtlasSizeTarget(16, 16);
  auto encodings = GenerateEncodings();
  auto atlaser = std::make_shared<FixedWidthAtlaser>(kAtlasSizeTarget);
  AtlasSizeTargetSelector encoding_selector(atlaser);
  std::vector<EncodedTexture> encoded_textures(encodings.size());
  encoding_selector.Select(encodings, absl::MakeSpan(encoded_textures));
  EXPECT_EQ(1, encoded_textures.size());
  EXPECT_GE(kAtlasSizeTarget[0], encoded_textures[0].image.Width());
  EXPECT_GE(kAtlasSizeTarget[1], encoded_textures[0].image.Height());
  EXPECT_NEAR(3456.0f, encoded_textures[0].rate, 1.0e-3f);
  EXPECT_NEAR(0.0f, encoded_textures[0].distortion, 0.1f);
}

}  // namespace

}  // namespace compressor
}  // namespace seurat
