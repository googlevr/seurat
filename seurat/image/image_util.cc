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

#include "seurat/image/image_util.h"

#include "ion/portgfx/glheaders.h"

namespace seurat {
namespace image {

namespace {

// This template function handles the actual conversion from Ion Image to
// Seurat Image.
template <int Format, int Type, typename PixelType>
Image<PixelType> DoConvert(const ion::gfx::ImagePtr& ion_image) {
  ion::gfx::Image::PixelFormat pixel_format =
      ion::gfx::Image::GetPixelFormat(ion_image->GetFormat());
  CHECK_EQ(Format, pixel_format.format);
  CHECK_EQ(Type, pixel_format.type);

  Image<PixelType> simage(ion_image->GetWidth(), ion_image->GetHeight());
  const auto data_container = ion_image->GetData();
  const PixelType* data = data_container->GetData<PixelType>();
  for (int y = 0; y < simage.Height(); ++y) {
    for (int x = 0; x < simage.Width(); ++x) {
      simage.At(x, y) = data[y * simage.Width() + x];
    }
  }
  return simage;
}

}  // namespace

ion::gfx::ImagePtr ConvertSeuratImageToIonImage(const Image4ui8& simage) {
  ion::gfx::ImagePtr img = base::CreateImage(ion::gfx::Image::Format::kRgba8888,
                                             {simage.Width(), simage.Height()});
  const auto data_container = img->GetData();
  base::Color4ui8* data = data_container->GetMutableData<base::Color4ui8>();
  for (int y = 0; y < simage.Height(); ++y) {
    for (int x = 0; x < simage.Width(); ++x) {
      data[y * simage.Width() + x] = simage.At(x, y);
    }
  }
  return img;
}

template <>
Image1f ConvertIonImageToSeuratImage<Image1f>(
    const ion::gfx::ImagePtr& ion_image) {
  return DoConvert<GL_RED, GL_FLOAT, base::Color1f>(ion_image);
}

template <>
Image4ui8 ConvertIonImageToSeuratImage<Image4ui8>(
    const ion::gfx::ImagePtr& ion_image) {
  return DoConvert<GL_RGBA, GL_UNSIGNED_BYTE, base::Color4ui8>(ion_image);
}

bool IsOpaque(const Image4f& image, float alpha_threshold) {
  for (int y = 0; y < image.Height(); ++y) {
    for (int x = 0; x < image.Width(); ++x) {
      if (image.At(x, y)[3] < alpha_threshold) return false;
    }
  }
  return true;
}

}  // namespace image
}  // namespace seurat
