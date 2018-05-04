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

#ifndef VR_SEURAT_COMPRESSOR_RGBA_RGBA_RATE_RESIZER_H_
#define VR_SEURAT_COMPRESSOR_RGBA_RGBA_RATE_RESIZER_H_

#include <memory>

#include "ion/math/vector.h"
#include "seurat/compressor/resampler/resampler.h"
#include "seurat/compressor/texture_encoder.h"

namespace seurat {
namespace compressor {

// RateResizer provides an implementation of the TextureEncoder interface.
// Textures are encoded by downsampling. For each downsampled version of each
// input texture the corresponding bitrate and distortion are computed and
// stored in an EncodedTexture data structure. Texture sizes are constrained
// to be multiples of a block size, matching the ASTC block size.
class RgbaRateResizer : public TextureEncoder {
 public:
  RgbaRateResizer(int thread_count, std::unique_ptr<Resampler> upsampler,
                  std::unique_ptr<Resampler> downsampler,
                  const ion::math::Vector2i& block_size);

  // Returns via |encodings| a vector of downsampled versions for each
  // texture in the input array |textures|.
  void Encode(absl::Span<const image::Image4f> textures,
              absl::Span<std::vector<EncodedTexture>> encodings) const override;

 private:
  // Number of threads.
  const int thread_count_;

  // Upsampler used to compute distortion.
  const std::unique_ptr<Resampler> upsampler_;

  // Downsampler used to compute the resized versions of textures.
  const std::unique_ptr<Resampler> downsampler_;

  // Block size used as texture size increments in x and y directions. Matches
  // the ASTC block size.
  const ion::math::Vector2i block_size_;
};

}  // namespace compressor
}  // namespace seurat

#endif  // VR_SEURAT_COMPRESSOR_RGBA_RGBA_RATE_RESIZER_H_
