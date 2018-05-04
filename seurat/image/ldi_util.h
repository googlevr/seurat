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

#ifndef VR_SEURAT_IMAGE_LDI_UTIL_H_
#define VR_SEURAT_IMAGE_LDI_UTIL_H_

#include "seurat/image/image.h"
#include "seurat/image/ldi.h"

namespace seurat {
namespace image {

// Flattens an |ldi| and returns the result as Image3f.  The samples in each
// pixel are composited using 'over' compositing. RGB values are assumed to
// be pre-multiplied by alpha.
Image3f FlattenLdi(const Ldi4f& ldi);

// Flips an |ldi| vertically and returns the result.
Ldi4f FlipLdiVertically(const Ldi4f& ldi);

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_LDI_UTIL_H_
