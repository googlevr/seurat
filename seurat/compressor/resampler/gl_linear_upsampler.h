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

#ifndef VR_SEURAT_COMPRESSOR_RESAMPLER_GL_LINEAR_UPSAMPLER_H_
#define VR_SEURAT_COMPRESSOR_RESAMPLER_GL_LINEAR_UPSAMPLER_H_

#include "ion/math/vector.h"
#include "seurat/compressor/resampler/resampler.h"
#include "seurat/image/image.h"

namespace seurat {
namespace compressor {

// Implements a method for upsampling images that simulates the OpenGL's
// GL_LINEAR algorithm. Half-pixel texture borders are assumed: the vertices of
// a textured quad have texture coordinates {0.5f, 0.5f},
// {tsize[0] - 0.5f, 0.5f}, {tsize[0] - 0.5f, tsize[1] - 0.5f},
// {0.5f, tsize[1] - 0.5f}, where tsize is the texture size.
class GlLinearUpsampler : public Resampler {
 public:
  // Returns an image of size |target_size| given the source image |image|.
  image::Image4f Resample(
      const image::Image4f& image,
      const ion::math::Vector2i& target_size) const override;
};

}  // namespace compressor
}  // namespace seurat

#endif  // VR_SEURAT_COMPRESSOR_RESAMPLER_GL_LINEAR_UPSAMPLER_H_
