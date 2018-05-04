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

#include <algorithm>
#include <functional>
#include <vector>

#include "ion/math/vector.h"
#include "seurat/base/parallel.h"
#include "seurat/base/util.h"
#include "seurat/compressor/image_metrics/image_metrics.h"
#include "seurat/compressor/resampler/resampler.h"
#include "seurat/image/image_util.h"

namespace seurat {
namespace compressor {

using image::Image4f;
using ion::math::Vector2i;

namespace {

// Returns a downsampled version of |image| of size |target_size| computed with
// the given |downsampler|. The bitrate |rate| and the image |distortion| are
// also returned.
Image4f ComputeResizedImage(const Image4f& image, const Resampler& upsampler,
                            const Resampler& downsampler,
                            const Vector2i& target_size, float* rate,
                            float* distortion) {
  Vector2i initial_size = image.GetSize();
  const ion::gfx::Image::Format kUncompressed =
      ion::gfx::Image::Format::kRgba8888;
  if (target_size == initial_size) {
    *rate = EvalBitrate(initial_size, kUncompressed);
    *distortion = 0.0f;
    return image;
  }
  Image4f resized_image = downsampler.Resample(image, target_size);
  Image4f distorted_image = upsampler.Resample(resized_image, initial_size);
  *rate = EvalBitrate(target_size, kUncompressed);
  *distortion = EvalSSD(image, distorted_image);
  return resized_image;
}

// Returns an array containing all resizings of |texture| allowed by the initial
// size of the texture and by the block size.
std::vector<EncodedTexture> ComputeAllResizings(const Image4f& texture,
                                                const Resampler& upsampler,
                                                const Resampler& downsampler,
                                                const Vector2i& block_size) {
  const Vector2i& initial_size = texture.GetSize();
  // The number of blocks that fit in each direction.
  Vector2i n_resizings((initial_size[0] + block_size[0] - 1) / block_size[0],
                       (initial_size[1] + block_size[1] - 1) / block_size[1]);
  // Minimum size of a downsampled texture. This constraint follows from the
  // half-texel border.
  constexpr int kMinimumTargetSize = 2;
  std::vector<EncodedTexture> all_resizings;
  for (int i = 0; i < n_resizings[0]; ++i) {
    for (int j = 0; j < n_resizings[1]; ++j) {
      Vector2i target_size((i + 1) * block_size[0], (j + 1) * block_size[1]);
      if (target_size[0] < kMinimumTargetSize ||
          target_size[1] < kMinimumTargetSize) {
        continue;
      }
      CHECK_LE(target_size[0], initial_size[0]);
      CHECK_LE(target_size[1], initial_size[1]);
      float rate;
      float distortion;
      Image4f resized_image = ComputeResizedImage(
          texture, upsampler, downsampler, target_size, &rate, &distortion);

      all_resizings.push_back({std::move(resized_image),
                               ion::gfx::Image::Format::kRgba8888, rate,
                               distortion});
    }
  }
  return all_resizings;
}

// Discards resized versions of a texture that do not reduce the distortion
// when increasing the bitrate.
std::vector<EncodedTexture> FilterResizings(
    std::vector<EncodedTexture>* raw_resizings) {
  // Order from low to high bitrates.
  std::sort(raw_resizings->begin(), raw_resizings->end(),
            [](const EncodedTexture& a, const EncodedTexture& b) {
              return a.rate < b.rate;
            });
  // Final list of resizing options.
  std::vector<EncodedTexture> resizings;
  // Minimum attained distortion
  float min_distortion = std::numeric_limits<float>::infinity();
  for (auto& resizing : *raw_resizings) {
    // Discard a resizing with a distortion not smaller than min_distortion.
    if (min_distortion <= resizing.distortion) {
      continue;
    }
    min_distortion = resizing.distortion;
    // Discard the last resizing in the array if the current resizing has the
    // same rate at a lower distortion.
    if (!resizings.empty() && resizing.rate == resizings.back().rate) {
      resizings.pop_back();
    }
    // Add the current resizing to the final list.
    resizings.push_back(std::move(resizing));
  }

  return resizings;
}

std::vector<EncodedTexture> ComputeSingleTextureResizingOptions(
    const Image4f& texture, const Resampler& upsampler,
    const Resampler& downsampler, const Vector2i& block_size) {
  std::vector<EncodedTexture> all_resizings =
      ComputeAllResizings(texture, upsampler, downsampler, block_size);
  return FilterResizings(&all_resizings);
}

}  // namespace

RgbaRateResizer::RgbaRateResizer(int thread_count,
                                 std::unique_ptr<Resampler> upsampler,
                                 std::unique_ptr<Resampler> downsampler,
                                 const Vector2i& block_size)
    : thread_count_(thread_count),
      upsampler_(std::move(upsampler)),
      downsampler_(std::move(downsampler)),
      block_size_(block_size) {
  CHECK_LT(0, block_size_[0]);
  CHECK_LT(0, block_size_[1]);
}

void RgbaRateResizer::Encode(
    absl::Span<const Image4f> textures,
    absl::Span<std::vector<EncodedTexture>> resizing_options) const {
  CHECK_EQ(textures.size(), resizing_options.size());
  int n_textures = textures.size();
  base::ParallelFor(thread_count_, n_textures, [&](int texture_index) {
    resizing_options[texture_index] = ComputeSingleTextureResizingOptions(
        textures[texture_index], *upsampler_, *downsampler_, block_size_);
  });
}

}  //  namespace compressor
}  //  namespace seurat
