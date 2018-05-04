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

#ifndef VR_SEURAT_API_LDI_EXPORTER_H_
#define VR_SEURAT_API_LDI_EXPORTER_H_

#include <string>
#include <vector>

#include "seurat/base/status.h"

namespace seurat {
namespace api {

// Writes a layered depth image (LDI) as Deep EXR to a string |file_contents|.
//
// An LDI is a image, where every pixel stores an arbitrary number of samples
// (can be zero). A sample consists of an RGBA color value and a depth
// value. Samples in each pixel must be sorted front to back. Pixels are laid
// out in scanline order with the first scanline at the bottom of the image.
//
// Example: Pixel (0,0) in the lower-left corner with samples a0, a1, a2; pixel
// (0,1) with samples b0, b1; pixel (1,0) with no samples, pixel (1, 1) with
// samples c0, c1.  Memory layout a0a1a2b0b1c0c1;
//
// |width| and |height| specify the dimensions of the LDI in pixels.
//
// |sample_counts| specifies the number of samples per pixel. The vector must
// have |width| x |height| elements laid out in scanline order with the first
// scanline at the bottom of the image.
//
// Samples are passed in two separate arrays |rgba_colors| and |depths| that
// follow the memory layout described above.
base::Status WriteDeepExr(int width, int height,
                          const std::vector<int>& sample_counts,
                          const std::vector<float>& rgba_colors,
                          const std::vector<float>& depths,
                          std::string* file_contents);

}  // namespace api
}  // namespace seurat

#endif  // VR_SEURAT_API_LDI_EXPORTER_H_
