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

#include "seurat/compressor/rgba/rgba_selecting_compressor.h"

#include <memory>
#include <numeric>
#include <vector>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/parallel.h"
#include "seurat/compressor/image_metrics/image_metrics.h"
#include "seurat/compressor/resampler/gl_linear_upsampler.h"
#include "seurat/compressor/rgba/rgba_compressor_util.h"
#include "seurat/image/fixed_width_atlaser.h"
#include "seurat/image/image.h"
#include "seurat/image/image_test_utils.h"
#include "seurat/image/image_util.h"

namespace seurat {
namespace compressor {

using base::Color4f;
using image::FixedWidthAtlaser;
using image::Image4f;
using ion::math::Point2i;
using ion::math::Vector2i;

namespace {

// Compute the total bitrate for the textures in the |textures| array.
float EvalTotalBitrate(const std::vector<Image4f>& textures) {
  return std::accumulate(textures.begin(), textures.end(), 0.0f,
                         [](float sum, const Image4f& texture) {
                           return sum + EvalBitrate(
                                            texture.GetSize(),
                                            ion::gfx::Image::Format::kRgba8888);
                         });
}

// Generates a set of checkerboard textures of the same size but increasing
// checker square sizes, so that higher index textures are easier to compress
// than lower index ones. Requests a maximum atlas size that will result in a 4x
// compression rate and verifies that the compressed textures have progressively
// lower sizes/bitrates.
TEST(RgbaSelectingCompressor, Compress) {
  const int kThreadCount = base::GetNumberOfHardwareThreads();
  const int kTextureCount = 4;
  const Vector2i kSourceSize = Vector2i(24, 32);
  const Vector2i kBlockSize(1, 1);
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kGreen(0.0f, 1.0f, 0.0f, 1.0f);
  const Vector2i kMaxAtlasSize(27, 40);
  const std::vector<Color4f> colors = {{kRed, kGreen}};

  std::vector<Image4f> textures;
  for (int i = 0; i < kTextureCount; ++i) {
    Image4f texture = image::MakeCheckerboard(kSourceSize, 4 * (i + 1), colors);
    textures.push_back(std::move(texture));
  }

  float total_uncompressed_rate = EvalTotalBitrate(textures);
  EXPECT_NEAR(98304.0f, total_uncompressed_rate, 1.0e-3f);

  auto atlaser = std::make_shared<FixedWidthAtlaser>(kMaxAtlasSize);
  auto rgba_compressor =
      BuildAtlasSizeTargetRgbaCompressor(kThreadCount, kBlockSize, atlaser);
  std::vector<Image4f> compressed_textures(kTextureCount);
  rgba_compressor->Compress(textures, absl::MakeSpan(compressed_textures));

  EXPECT_NEAR(27648.0f, EvalTotalBitrate(compressed_textures), 1.0e-3f);
  std::vector<float> expected_rates = {13824.0f, 6144.0f, 4608.0f, 3072.0f};
  std::vector<Vector2i> expected_sizes = {{18, 24}, {12, 16}, {6, 24}, {12, 8}};

  // Check the size of the atlas obtained by laying out the compressed textures.
  std::vector<Vector2i> texture_sizes;
  for (const auto& texture : compressed_textures) {
    texture_sizes.push_back(texture.GetSize());
  }
  Vector2i atlas_size;
  std::vector<Point2i> layout(kTextureCount);
  atlaser->LayoutTiles(texture_sizes, &atlas_size, absl::MakeSpan(layout));
  EXPECT_GE(kMaxAtlasSize[0], atlas_size[0]);
  EXPECT_GE(kMaxAtlasSize[1], atlas_size[1]);

  for (int i = 0; i < kTextureCount; ++i) {
    EXPECT_NEAR(expected_rates[i],
                EvalBitrate(compressed_textures[i].GetSize(),
                            ion::gfx::Image::Format::kRgba8888),
                1.0e-3f)
        << " i " << i;
    EXPECT_EQ(expected_sizes[i], compressed_textures[i].GetSize())
        << " i " << i;
  }
}

}  // namespace

}  // namespace compressor
}  // namespace seurat
