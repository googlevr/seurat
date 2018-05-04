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

#include "seurat/image/rgbaui8_codec.h"

#include "seurat/base/array2d_util.h"
#include "seurat/image/image_util.h"

namespace seurat {
namespace image {

ion::gfx::ImagePtr RgbaUi8Codec::Compress(const Image4f& simage) const {
  return ConvertSeuratImageToIonImage(simage);
}

Image4f RgbaUi8Codec::Decompress(const ion::gfx::ImagePtr& ion_image) const {
  Image4ui8 simage = ConvertIonImageToSeuratImage<Image4ui8>(ion_image);
  Image4f simage4f(simage.Width(), simage.Height());
  base::TransformArray(simage, &simage4f,
                       [](const base::Color4ui8& p) { return p.AsColorF(); });
  return simage4f;
}

}  // namespace image
}  // namespace seurat
