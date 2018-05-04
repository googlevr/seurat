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

#ifndef VR_SEURAT_COMPRESSOR_RGBA_RGBA_SELECTING_COMPRESSOR_H_
#define VR_SEURAT_COMPRESSOR_RGBA_RGBA_SELECTING_COMPRESSOR_H_

#include <memory>

#include "seurat/compressor/rgba/rgba_compressor.h"
#include "seurat/compressor/rgba/rgba_encoding_selector.h"
#include "seurat/compressor/texture_encoder.h"
#include "seurat/image/image.h"

namespace seurat {
namespace compressor {

// RgbaSelectorCompressor specifies an implementation of the RgbaCompressor
// interface based on making selections from sets of texture encodings.
class RgbaSelectingCompressor : public RgbaCompressor {
 public:
  RgbaSelectingCompressor(
      std::unique_ptr<TextureEncoder> texture_encoder,
      std::unique_ptr<RgbaEncodingSelector> encoding_selector);

  // Returns via |compressed_textures| one downsampled version for each of the
  // input |textures|.
  void Compress(absl::Span<const image::Image4f> textures,
                absl::Span<image::Image4f> compressed_textures) const override;

 private:
  // Generates several encodings for each of the input textures.
  const std::unique_ptr<TextureEncoder> texture_encoder_;

  // Selects one encoding for each input texture out of the encodings built by
  // the texture encoder.
  const std::unique_ptr<RgbaEncodingSelector> encoding_selector_;
};

}  // namespace compressor
}  // namespace seurat

#endif  // VR_SEURAT_COMPRESSOR_RGBA_RGBA_SELECTING_COMPRESSOR_H_
