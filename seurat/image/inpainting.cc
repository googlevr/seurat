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

#include "seurat/image/inpainting.h"

#include <array>
#include <vector>

#include "ion/math/vector.h"

namespace seurat {
namespace image {

using base::Array2D;
using base::Color4f;
using ion::math::Point2i;
using ion::math::Vector2i;

namespace {

// Inpaints the masked image by performing several iterations.
void InpaintSmoothSingleResolution(const Array2D<bool>& mask, Image4f* image,
                                   int iteration_count) {
  // All masked pixel coordinates.
  std::vector<Point2i> masked;
  for (int y = 0; y < mask.Height(); ++y) {
    for (int x = 0; x < mask.Width(); ++x) {
      if (mask.At(x, y)) {
        masked.emplace_back(x, y);
      }
    }
  }

  // 4-connectivity neighborhood.
  const std::array<Vector2i, 4> neighborhood = {{
      {0, -1},  //
      {-1, 0},  //
      {1, 0},   //
      {0, 1}    //
  }};
  // Perform several iterations, looping over all pixels to be solved for.
  //
  // This is equivalent to using Gauss-Seidel to solve the linear system
  // corresponding to finding the image which minimizes its laplacian.
  //
  // Each iteration simply replaces each pixel with the average of its
  // neighbors.
  for (int i = 0; i < iteration_count; ++i) {
    for (const auto& coords : masked) {
      int neighbors = 0;
      Color4f neighbor_sum = Color4f::Zero();
      for (const auto& offset : neighborhood) {
        const Point2i neighbor_coords = coords + offset;
        if (!image->IsInside(neighbor_coords)) {
          continue;
        }
        neighbor_sum += image->At(neighbor_coords);
        neighbors++;
      }
      if (neighbors == 0) continue;

      image->At(coords[0], coords[1]) =
          neighbor_sum / static_cast<float>(neighbors);
    }
  }
}

// Downsample using a 2x2 box filter, with the filter only evaluated at
// *unmasked* samples.
void MaskedDownsample(const Image4f& original, const Array2D<bool>& mask,
                      Image4f* downsampled) {
  for (int y = 0; y < downsampled->Height(); ++y) {
    for (int x = 0; x < downsampled->Width(); ++x) {
      downsampled->At(x, y) = Color4f::Zero();
      int unmasked_sample_count = 0;
      for (int v = 0; v < 2; ++v) {
        for (int u = 0; u < 2; ++u) {
          Point2i coord(x * 2 + u, y * 2 + v);
          if (coord[0] >= mask.Width()) {
            coord[0] = mask.Width() - 1;
          }
          if (coord[1] >= mask.Height()) {
            coord[1] = mask.Height() - 1;
          }
          if (!mask.At(coord)) {
            downsampled->At(x, y) += original.At(coord);
            unmasked_sample_count++;
          }
        }
      }
      if (unmasked_sample_count > 0) {
        downsampled->At(x, y) /= unmasked_sample_count;
      }
    }
  }
}

// Overwrites *masked* pixels of |upsampled| with the |original| image.
void MaskedUpsample(const Image4f& original, const Array2D<bool>& mask,
                    Image4f* upsampled) {
  CHECK_EQ(upsampled->GetSize(), mask.GetSize());
  CHECK_EQ((upsampled->GetSize() + Vector2i(1, 1)) / 2, original.GetSize());
  for (int y = 0; y < upsampled->Height(); ++y) {
    for (int x = 0; x < upsampled->Width(); ++x) {
      if (mask.At(x, y)) {
        upsampled->At(x, y) = original.At(x / 2, y / 2);
      }
    }
  }
}

// Downsample using the logical-and of a 2x2 neighborhood.
void Downsample(const Array2D<bool>& original, Array2D<bool>* downsampled) {
  for (int y = 0; y < downsampled->Height(); ++y) {
    for (int x = 0; x < downsampled->Width(); ++x) {
      downsampled->At(x, y) = true;
      for (int v = 0; v < 2; ++v) {
        for (int u = 0; u < 2; ++u) {
          if (!original.IsInside({x * 2 + u, y * 2 + v})) {
            continue;
          }
          downsampled->At(x, y) &= original.At(x * 2 + u, y * 2 + v);
        }
      }
    }
  }
}

}  // namespace

void InpaintSmooth(const base::Array2D<bool>& mask, image::Image4f* image) {
  // Operate in a coarse-to-fine manner over a pyramid, performing several
  // iterations at each level.

  const int kIterationCount = 30;

  CHECK_EQ(image->GetSize(), mask.GetSize());
  if (image->GetSize()[0] > 2 && image->GetSize()[1] > 2) {
    Image4f downsampled_rgba((image->GetSize() + Vector2i(1, 1)) / 2);
    MaskedDownsample(*image, mask, &downsampled_rgba);

    Array2D<bool> downsampled_mask((image->GetSize() + Vector2i(1, 1)) / 2);
    Downsample(mask, &downsampled_mask);

    // Recurse on downsampled problem.
    InpaintSmooth(downsampled_mask, &downsampled_rgba);

    // Upsample the result.
    MaskedUpsample(downsampled_rgba, mask, image);
  }

  InpaintSmoothSingleResolution(mask, image, kIterationCount);
}

}  // namespace image
}  // namespace seurat
