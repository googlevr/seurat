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

#include "seurat/viewer/viewer_util.h"

#include "ion/gfx/framebufferobject.h"
#include "ion/gfx/image.h"
#include "ion/image/conversionutils.h"

namespace seurat {
namespace viewer {

std::string CaptureFrame(const ion::gfx::RendererPtr& renderer) {
  const ion::gfx::FramebufferObjectPtr framebuffer =
      renderer->GetCurrentFramebuffer();

  const ion::gfx::ImagePtr image = renderer->ReadImage(
      ion::math::Range2i(ion::math::Point2i(0, 0),
                         ion::math::Point2i(framebuffer->GetWidth(),
                                            framebuffer->GetHeight())),
      ion::gfx::Image::kRgb888, ion::base::AllocatorPtr());

  // OpenGL / ion and image / PNG coordinates are different. This is
  // necessary for the saved image to be like what's rendered on screen.
  const bool flip_vertically = true;
  const std::vector<uint8> png_image = ion::image::ConvertToExternalImageData(
      image, ion::image::kPng, flip_vertically);

  return std::string(reinterpret_cast<const char*>(png_image.data()),
                     png_image.size());
}

}  // namespace viewer
}  // namespace seurat
