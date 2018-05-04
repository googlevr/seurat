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

#ifndef VR_SEURAT_IMAGE_RGBAUI8_CODEC_H_
#define VR_SEURAT_IMAGE_RGBAUI8_CODEC_H_

#include "ion/gfx/image.h"
#include "ion/math/vector.h"
#include "seurat/image/codec.h"
#include "seurat/image/image.h"

namespace seurat {
namespace image {

// RgbaUi8Codec compresses Image4f by quantizing to RGBA8 format.
class RgbaUi8Codec : public Codec {
 public:
  ~RgbaUi8Codec() override = default;

  ion::gfx::Image::Format GetFormat() const override {
    return ion::gfx::Image::Format::kRgba8888;
  }

  ion::math::Vector2i GetBlockSize() const override {
    return ion::math::Vector2i::Zero();
  }
  ion::gfx::ImagePtr Compress(const image::Image4f& simage) const override;
  image::Image4f Decompress(const ion::gfx::ImagePtr& ion_image) const override;
};

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_RGBAUI8_CODEC_H_
