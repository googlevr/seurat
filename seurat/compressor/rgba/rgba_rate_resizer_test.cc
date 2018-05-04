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

#include "seurat/compressor/rgba/rgba_rate_resizer.h"

#include <vector>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/parallel.h"
#include "seurat/base/util.h"
#include "seurat/compressor/image_metrics/image_metrics.h"
#include "seurat/compressor/resampler/box_downsampler.h"
#include "seurat/compressor/resampler/gl_linear_upsampler.h"
#include "seurat/image/image.h"
#include "seurat/image/image_test_utils.h"
#include "seurat/image/image_util.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace compressor {

using base::Color4f;
using image::Image4f;
using ion::math::Vector2i;

namespace {

// Compute all resizing options for a checkerboard texture and check the number
// and sizes or downsampled images.
TEST(RgbaRateResizer, Encode) {
  constexpr int kSquareSize = 2;
  const int kThreadCount = base::GetNumberOfHardwareThreads();
  const Vector2i kSourceSize = Vector2i(9, 12);
  const Vector2i kBlockSize(1, 1);
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kGreen(0.0f, 1.0f, 0.0f, 1.0f);
  const Color4f kBlue(0.0f, 0.0f, 1.0f, 1.0f);
  const std::vector<Color4f> colors = {{kRed, kGreen, kBlue}};
  Image4f texture = image::MakeCheckerboard(kSourceSize, kSquareSize, colors);
  RgbaRateResizer resizer(
      kThreadCount, std::unique_ptr<Resampler>(new GlLinearUpsampler()),
      std::unique_ptr<Resampler>(new BoxDownsampler()), kBlockSize);
  std::vector<std::vector<EncodedTexture>> resizings(1);
  resizer.Encode({texture}, absl::MakeSpan(resizings));
  int option_count = resizings[0].size();
  std::vector<Vector2i> expected_sizes = {
      {{2, 2}, {3, 2},  {4, 2},  {5, 2},  {3, 4},  {3, 5}, {4, 4}, {3, 6},
       {4, 5}, {4, 6},  {5, 6},  {6, 6},  {7, 6},  {8, 6}, {9, 6}, {6, 12},
       {9, 9}, {7, 12}, {9, 10}, {8, 12}, {9, 11}, {9, 12}}};
  EXPECT_EQ(expected_sizes.size(), option_count);

  for (int j = 0; j < option_count; ++j) {
    EXPECT_EQ(expected_sizes[j], resizings[0][j].image.GetSize());
    if (j > 0) {
      EXPECT_LE(resizings[0][j - 1].rate, resizings[0][j].rate);
      EXPECT_LT(resizings[0][j].distortion, resizings[0][j - 1].distortion);
    }
  }
}

}  // namespace

}  // namespace compressor
}  // namespace seurat
