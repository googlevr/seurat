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

#ifndef VR_SEURAT_IMAGE_INPAINTING_H_
#define VR_SEURAT_IMAGE_INPAINTING_H_

#include "seurat/base/array2d.h"
#include "seurat/image/image.h"

namespace seurat {
namespace image {

// Inpaints the given |image| wherever |mask| is set to 'true'.
//
// Notably, the goal is *not* to generate a plausible/realistic inpainted
// texture.  Instead, the goal is to inpaint while maximizing smoothness
// (minimizing the laplacian of the resulting image).
//
// Works best for square images.
void InpaintSmooth(const base::Array2D<bool>& mask, image::Image4f* image);

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_INPAINTING_H_
