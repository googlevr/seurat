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

#ifndef VR_SEURAT_IMAGE_CODEC_H_
#define VR_SEURAT_IMAGE_CODEC_H_

#include "ion/gfx/image.h"
#include "ion/math/vector.h"
#include "seurat/image/image.h"

namespace seurat {
namespace image {

// Codec defines the interface for compression from Seurat Image4f to an Ion
// texture format, and decompression of the reverse (say, for analysis of
// compression artifacts).
class Codec {
 public:
  virtual ~Codec() = default;

  // Returns the Ion texture format.
  virtual ion::gfx::Image::Format GetFormat() const = 0;

  // Retrieves the pixel block granularity the compression imposes on input
  // images. Some (Ion) Image formats compress in blocks, e.g. DXT requires 4x4
  // blocks, and so require image dimensions be a multiple of this size.
  // Compressors not requiring such block granularity return zero for both
  // dimensions.
  virtual ion::math::Vector2i GetBlockSize() const = 0;

  // Compresses a Seurat |image| into an Ion image.
  virtual ion::gfx::ImagePtr Compress(const image::Image4f& image) const = 0;

  // Decompresses an Ion |image| into a Seurat image. The Ion image must be
  // in a format produced by the Compress API.
  virtual image::Image4f Decompress(const ion::gfx::ImagePtr& image) const = 0;
};

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_CODEC_H_
