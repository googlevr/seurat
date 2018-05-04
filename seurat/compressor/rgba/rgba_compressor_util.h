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

#ifndef VR_SEURAT_COMPRESSOR_RGBA_RGBA_COMPRESSOR_UTIL_H_
#define VR_SEURAT_COMPRESSOR_RGBA_RGBA_COMPRESSOR_UTIL_H_

#include <memory>

#include "ion/math/vector.h"
#include "seurat/compressor/rgba/rgba_compressor.h"
#include "seurat/image/atlaser.h"

namespace seurat {
namespace compressor {

// Builds an RgbaCompressor instance that will minimize the total distortion
// of textures given a texture atlasing process specified by |atlaser|. The
// atlasing process may include constraints on the size of the resulting
// texture atlas.
std::unique_ptr<RgbaCompressor> BuildAtlasSizeTargetRgbaCompressor(
    int thread_count, const ion::math::Vector2i& texture_alignment,
    std::shared_ptr<image::Atlaser> atlaser);

}  // namespace compressor
}  // namespace seurat

#endif  // VR_SEURAT_COMPRESSOR_RGBA_RGBA_COMPRESSOR_UTIL_H_
