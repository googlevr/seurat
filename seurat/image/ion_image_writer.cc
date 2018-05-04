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

#include "seurat/image/ion_image_writer.h"

#include "ion/image/conversionutils.h"
#include "absl/strings/substitute.h"
#include "seurat/image/image_util.h"

namespace seurat {
namespace image {

base::Status WriteImagePngToString(const Image4f& image, std::string* encoded) {
  if (image.Width() <= 0 || image.Height() <= 0) {
    return base::InvalidArgumentError(absl::Substitute(
        "Image size ($0, $1) is not positive", image.Width(), image.Height()));
  }
  ion::gfx::ImagePtr ion_image = ConvertSeuratImageToIonImage(image);
  // Don't flip. Image flipping is handled by exporters.
  const bool kIonFlipImageVertically = false;
  std::vector<uint8> encoded_png = ion::image::ConvertToExternalImageData(
      ion_image, ion::image::ExternalImageFormat::kPng,
      kIonFlipImageVertically);
  encoded->resize(encoded_png.size());
  std::copy(encoded_png.begin(), encoded_png.end(), encoded->begin());
  return base::OkStatus();
}

}  // namespace image
}  // namespace seurat
