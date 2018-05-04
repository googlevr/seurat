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

#ifndef VR_SEURAT_IMAGE_IMAGE_UTIL_H_
#define VR_SEURAT_IMAGE_IMAGE_UTIL_H_

#include <utility>

#include "ion/gfx/image.h"
#include "seurat/base/color.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/image/image.h"

namespace seurat {
namespace image {

// Maps a seurat image type to the corresponding 8-bit integer pixel type and
// ion image format.
template <typename ImageType>
struct ConvertImageTraits {};

template <>
struct ConvertImageTraits<image::Image1f> {
  using IntPixelType = base::Color1ui8;
  static constexpr ion::gfx::Image::Format kIonFormat =
      ion::gfx::Image::Format::kAlpha;
};

template <>
struct ConvertImageTraits<image::Image3f> {
  using IntPixelType = base::Color3ui8;
  static constexpr ion::gfx::Image::Format kIonFormat =
      ion::gfx::Image::Format::kRgb888;
};

template <>
struct ConvertImageTraits<image::Image4f> {
  using IntPixelType = base::Color4ui8;
  static constexpr ion::gfx::Image::Format kIonFormat =
      ion::gfx::Image::Format::kRgba8888;
};

// Converts a floating-point Image to the corresponding 8 bits per channel ion
// image format.
template <typename ImageType>
ion::gfx::ImagePtr ConvertSeuratImageToIonImage(const ImageType& simage) {
  using IntPixelType = typename ConvertImageTraits<ImageType>::IntPixelType;
  ion::gfx::ImagePtr img =
      base::CreateImage(ConvertImageTraits<ImageType>::kIonFormat,
                        {simage.Width(), simage.Height()});
  const auto data_container = img->GetData();
  IntPixelType* data = data_container->GetMutableData<IntPixelType>();
  for (int y = 0; y < simage.Height(); ++y) {
    for (int x = 0; x < simage.Width(); ++x) {
      data[y * simage.Width() + x] = simage.At(x, y).AsColorUI8();
    }
  }
  return img;
}

// Converts an 8-bit Image to an ion image in rgba8888 format.
ion::gfx::ImagePtr ConvertSeuratImageToIonImage(const Image4ui8& simage);

// Converts an Ion Image to a Seurat Image.  This templated function is only
// specialized for the supported destination formats.
template <typename ReturnImageType>
ReturnImageType ConvertIonImageToSeuratImage(
    const ion::gfx::ImagePtr& ion_image);

// Returns true if |image| is opaque relative to the |alpha_threshold|.
bool IsOpaque(const Image4f& image, float alpha_threshold);

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_IMAGE_UTIL_H_
