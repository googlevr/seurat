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

#include "seurat/compressor/resampler/gl_linear_upsampler.h"

#include "ion/math/matrix.h"
#include "ion/math/range.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "seurat/image/image.h"

namespace seurat {
namespace compressor {

using image::Image;
using image::Image4f;
using ion::math::Matrix3f;
using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Range2i;
using ion::math::Vector2f;
using ion::math::Vector2i;

namespace {

// Returns the matrix of the source texture space from target texture space
// transform.
Matrix3f SourceFromTargetTexture(const Vector2i& source_size,
                                 const Vector2i& target_size) {
  Vector2f scale((static_cast<float>(source_size[0]) - 1.0f) /
                     (static_cast<float>(target_size[0]) - 1.0f),
                 (static_cast<float>(source_size[1]) - 1.0f) /
                     (static_cast<float>(target_size[1]) - 1.0f));
  return Matrix3f(scale[0], 0.0f, 0.5f * (1.0f - scale[0]),  //
                  0.0f, scale[1], 0.5f * (1.0f - scale[1]),  //
                  0.0f, 0.0f, 1.0f);
}

// Returns the barycentric coordinates for the sample point defined by
// |texcoords|. Simulating OpenGL's GL_LINEAR behavior, the barycentric
// coordinates are relative to a 2*2 block of texels in which the lower-left
// texel is indexed by the returned value |texel_index|. The texture coordinates
// |texcoords| span a texture of size |texture_size|.
Point2f ComputeSampleBary(const Point2f& texcoords,
                          const Vector2i& texture_size, Point2i* texel_index) {
  CHECK_LE(0.0f, texcoords[0]);
  CHECK_LE(texcoords[0], static_cast<float>(texture_size[0]));
  CHECK_LE(0.0f, texcoords[1]);
  CHECK_LE(texcoords[1], static_cast<float>(texture_size[1]));
  CHECK_LE(2, texture_size[0]);
  CHECK_LE(2, texture_size[1]);

  // Texel's x and y indices, in the range [0, texture_size-1].
  (*texel_index)[0] = static_cast<int>(std::floor(texcoords[0] - 0.5f));
  (*texel_index)[1] = static_cast<int>(std::floor(texcoords[1] - 0.5f));

  // |texel_index| must be in the interval [0, texture_size-2].
  (*texel_index)[0] =
      std::max(0, std::min(texture_size[0] - 2, (*texel_index)[0]));
  (*texel_index)[1] =
      std::max(0, std::min(texture_size[1] - 1, (*texel_index)[1]));

  Point2f bary = {texcoords[0] - 0.5f - (*texel_index)[0],
                  texcoords[1] - 0.5f - (*texel_index)[1]};

  for (int i = 0; i < 2; ++i) {
    if (bary[i] < 0.0f) {
      bary[i] = 0.0f;
    }
    if (1.0f < bary[i]) {
      bary[i] = 1.0f;
    }
  }

  return bary;
}

// Adds to the |upsampled_image| at position |upsampled_texel_index| the
// contribution of the source image texel with index |texel_index| and weight
// |weight|.
template <typename T>
void AddInTexel(const Image<T>& image, const Point2i& texel_index, float weight,
                const Point2i& upsampled_texel_index,
                Image<T>* upsampled_image) {
  Range2i image_range(Point2i(0, 0),
                      Point2i(0, 0) + image.GetSize() - Vector2i(1, 1));
  if (image_range.ContainsPoint(texel_index)) {
    upsampled_image->At(upsampled_texel_index) +=
        weight * image.At(texel_index);
  }
}

// Returns an upsampled version of |image| computed using a procedure that
// simulates the OpenGL's GL_LINEAR algorithm. Half-pixel texture borders are
// assumed: the vertices of a textured quad have texture coordinates
// {{0.5f, 0.5f}, {tsize[0] - 0.5f, 0.5f}, {tsize[0] - 0.5f, tsize[1] - 0.5f},
// {0.5f, tsize[1] - 0.5f}}, where tsize is the texture size.
template <typename T>
Image<T> Upsample(const Image<T>& image, const Vector2i& target_size) {
  Vector2i image_size = image.GetSize();
  CHECK_NE(image_size, target_size);
  CHECK_LE(image.Width(), target_size[0]);
  CHECK_LE(image.Height(), target_size[1]);
  Image<T> upsampled_image(target_size, T::Zero());
  for (int y = 0; y < target_size[1]; ++y) {
    for (int x = 0; x < target_size[0]; ++x) {
      Point2i upsampled_texel_index(x, y);
      // Compute the source texture coordinates of sample at (x,y) in target
      // texture space.
      Matrix3f source_from_target_texture =
          SourceFromTargetTexture(image_size, target_size);
      Point2f target_texcoords(x + 0.5f, y + 0.5f);
      Point2f texcoords = source_from_target_texture * target_texcoords;
      Point2i texel_index;
      Point2f bary = ComputeSampleBary(texcoords, image_size, &texel_index);

      // Add in contributions from texels belonging to the 2*2 block of texels
      // in the GL_LINEAR algorithm.

      // Lower left texel.
      AddInTexel(image, texel_index, (1.0f - bary[0]) * (1.0f - bary[1]),
                 upsampled_texel_index, &upsampled_image);
      // Lower right texel.
      AddInTexel(image, texel_index + Vector2i(1, 0),
                 bary[0] * (1.0f - bary[1]), upsampled_texel_index,
                 &upsampled_image);
      // Upper left texel.
      AddInTexel(image, texel_index + Vector2i(0, 1),
                 (1.0f - bary[0]) * bary[1], upsampled_texel_index,
                 &upsampled_image);
      // Upper right texel.
      AddInTexel(image, texel_index + Vector2i(1, 1), bary[0] * bary[1],
                 upsampled_texel_index, &upsampled_image);

      for (int c = 0; c < T::kDimension; ++c) {
        upsampled_image.At(x, y)[c] =
            std::max(0.0f, std::min(1.0f, upsampled_image.At(x, y)[c]));
      }
    }
  }
  return upsampled_image;
}

}  // namespace

Image4f GlLinearUpsampler::Resample(const Image4f& image,
                                    const Vector2i& target_size) const {
  return Upsample(image, target_size);
}

}  // namespace compressor
}  // namespace seurat
