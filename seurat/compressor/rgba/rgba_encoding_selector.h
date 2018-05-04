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

#ifndef VR_SEURAT_COMPRESSOR_RGBA_RGBA_ENCODING_SELECTOR_H_
#define VR_SEURAT_COMPRESSOR_RGBA_RGBA_ENCODING_SELECTOR_H_

#include <memory>
#include <vector>

#include "absl/types/span.h"
#include "seurat/compressor/texture_encoder.h"
#include "seurat/image/atlaser.h"
#include "seurat/image/image.h"

namespace seurat {
namespace compressor {

// Interface for selecting texture encodings.
class RgbaEncodingSelector {
 public:
  virtual ~RgbaEncodingSelector() = default;

  // Returns via the |encoded_textures| array one encoding for each texture
  // which was selected from the set of |encodings|. Both the |encodings| and
  // |encoded_textures| arrays are of the same size, which equals the number of
  // textures. Each entry in the |encodings| array is itself an array containing
  // all texture encodings that belong to the same texture.
  virtual void Select(absl::Span<const std::vector<EncodedTexture>> encodings,
                      absl::Span<EncodedTexture> encoded_textures) = 0;
};

// Implementation of the RgbaEncodingSelector interface based on Lagrange
// optimization. Distortion due to texture encoding is minimized subject to a
// constraint on the size of the texture atlas generated from the selected
// textures.
class AtlasSizeTargetSelector : public RgbaEncodingSelector {
 public:
  explicit AtlasSizeTargetSelector(std::shared_ptr<image::Atlaser> atlaser);

  void Select(absl::Span<const std::vector<EncodedTexture>> encodings,
              absl::Span<EncodedTexture> encoded_textures) override;

 private:
  // Defines the process of laying out texture tiles in a texture atlas,
  // including possible constraints on the atlas size.
  std::shared_ptr<image::Atlaser> atlaser_;
};

}  // namespace compressor
}  // namespace seurat

#endif  // VR_SEURAT_COMPRESSOR_RGBA_RGBA_ENCODING_SELECTOR_H_
