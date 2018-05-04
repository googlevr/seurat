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
#include <vector>

#include "seurat/base/util.h"
#include "seurat/compressor/rgba/rgba_compressor.h"
#include "seurat/compressor/rgba/rgba_encoding_selector.h"
#include "seurat/compressor/texture_encoder.h"
#include "seurat/image/image.h"

namespace seurat {
namespace compressor {

using image::Image4f;

RgbaSelectingCompressor::RgbaSelectingCompressor(
    std::unique_ptr<TextureEncoder> texture_encoder,
    std::unique_ptr<RgbaEncodingSelector> encoding_selector)
    : texture_encoder_(std::move(texture_encoder)),
      encoding_selector_(std::move(encoding_selector)) {}

void RgbaSelectingCompressor::Compress(
    absl::Span<const Image4f> textures,
    absl::Span<Image4f> compressed_textures) const {
  CHECK_EQ(textures.size(), compressed_textures.size());
  int n_textures = textures.size();
  std::vector<std::vector<EncodedTexture>> encoding_options(n_textures);
  texture_encoder_->Encode(textures, absl::MakeSpan(encoding_options));
  std::vector<EncodedTexture> encoded_textures(n_textures);
  encoding_selector_->Select(encoding_options,
                             absl::MakeSpan(encoded_textures));
  for (int i = 0; i < n_textures; ++i) {
    compressed_textures[i] = std::move(encoded_textures[i].image);
  }
}

}  // namespace compressor
}  // namespace seurat
