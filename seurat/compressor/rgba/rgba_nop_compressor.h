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

#ifndef VR_SEURAT_COMPRESSOR_RGBA_RGBA_NOP_COMPRESSOR_H_
#define VR_SEURAT_COMPRESSOR_RGBA_RGBA_NOP_COMPRESSOR_H_

#include "seurat/compressor/rgba/rgba_compressor.h"

namespace seurat {
namespace compressor {

// RgbaNopCompressor implements the RgbaCompressor interface by just copying the
// input textures to output.
class RgbaNopCompressor : public RgbaCompressor {
 public:
  // Returns via |compressed_textures| one copy for each of the input
  // |textures|.
  void Compress(absl::Span<const image::Image4f> textures,
                absl::Span<image::Image4f> compressed_textures) const override;
};

}  // namespace compressor
}  // namespace seurat

#endif  // VR_SEURAT_COMPRESSOR_RGBA_RGBA_NOP_COMPRESSOR_H_
