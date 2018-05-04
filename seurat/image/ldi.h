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

#ifndef VR_SEURAT_IMAGE_LDI_H_
#define VR_SEURAT_IMAGE_LDI_H_

#include <numeric>
#include <vector>

#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/base/array2d.h"
#include "seurat/base/color.h"

namespace seurat {
namespace image {

// A class template for layered depth images (LDIs). Only the color values of
// samples are mutable. Size, number of samples in each pixel, and the depth
// values are immutable. This is required for an efficient, dense memory layout.
//
// Color and depth is laid out in separate, dense arrays in memory. Samples of
// each pixel are laid out contiguously in memory. Pixels are laid out in
// scanline order. Example: Pixel (0,0) with samples a0, a1, a2; pixel (0,1)
// with samples b0, b1; pixel (1,0) with no samples, pixel (1, 1) with samples
// c0, c1.  Memory layout a0a1a2b0b1c0c1;
//
// Samples in each pixel (color and depth) are sorted front to back according to
// the depth samples.
template <typename ColorType>
class Ldi {
 public:
  Ldi() : offsets_{0} {}

  // Construct from deep image data with an arbitrary number of samples per
  // pixel (including zero).
  Ldi(const ion::math::Vector2i& size, std::vector<int> sample_counts,
      std::vector<ColorType> colors, std::vector<float> depths);

  // Construct from image data with exactly one sample per pixel.
  Ldi(const ion::math::Vector2i& size, std::vector<ColorType> colors,
      std::vector<float> depths);

  // Returns the dimensions of the LDI in pixels.
  const ion::math::Vector2i& GetSize() const { return size_; }

  // Returns the width of the LDI in pixels.
  int GetWidth() const { return size_[0]; }

  // Returns the height of the LDI in pixels.
  int GetHeight() const { return size_[1]; }

  // Returns the number of samples in the given |pixel|.
  int GetSampleCount(const ion::math::Point2i& pixel) const;

  // Return the number of samples in all pixels.
  int GetSampleCount() const;

  // Returns a span for the color samples in the given |pixel|.
  absl::Span<const ColorType> GetColors(const ion::math::Point2i& pixel) const;

  // Returns a span for the color samples in the given |pixel|.
  absl::Span<ColorType> GetMutableColors(const ion::math::Point2i& pixel);

  // Returns a span for the depth samples in the given |pixel|.
  absl::Span<const float> GetDepths(const ion::math::Point2i& pixel) const;

  // Returns a span for all color values in the LDI.
  absl::Span<const ColorType> GetColors() const {
    return absl::Span<const ColorType>(colors_.data(), colors_.size());
  }

  // Returns a span for all color values in the LDI.
  absl::Span<ColorType> GetMutableColors() {
    return absl::Span<ColorType>(colors_.data(), colors_.size());
  }

  // Returns a span for all depth values in the LDI.
  absl::Span<const float> GetDepths() const {
    return absl::Span<const float>(depths_.data(), depths_.size());
  }

  // Computes the offset and the |sample_count| for a given |pixel|.
  void GetSampleCountAndOffset(const ion::math::Point2i& pixel, int* offset,
                               int* sample_count) const;

 private:
  // Size of the LDI in pixels.
  ion::math::Vector2i size_;

  // A dense, flat array storing the color values for all pixels. See class
  // comment for the memory layout.
  std::vector<ColorType> colors_;

  // A dense, flat array storing the depth values for all pixels. See class
  // comment for the memory layout.
  std::vector<float> depths_;

  // Offsets into the color and depth vectors for each pixel, laid out in
  // scanline order.
  //
  // The vector is terminated by a sentinel, such that the sample count can
  // always be computed as the difference to the next offset.
  std::vector<int> offsets_;
};

template <typename ColorType>
Ldi<ColorType>::Ldi(const ion::math::Vector2i& size,
                    std::vector<int> sample_counts,
                    std::vector<ColorType> colors, std::vector<float> depths)
    : size_(size),
      colors_(std::move(colors)),
      depths_(std::move(depths)),
      offsets_(std::move(sample_counts)) {
  CHECK_GE(size_[0], 0);
  CHECK_GE(size_[1], 0);
  CHECK_EQ(colors_.size(), depths_.size());
  CHECK_EQ(offsets_.size(), size_[0] * size_[1]);
  // Reserve space for the sentinel. Ideally the caller has done that already.
  offsets_.reserve(size_[0] * size_[1] + 1);
  // Compute the offsets as the prefix sum of the sample counts in place.
  int current_offset = 0;
  for (int i = 0; i < offsets_.size(); ++i) {
    int current_sample_count = offsets_[i];
    offsets_[i] = current_offset;
    current_offset += current_sample_count;
  }
  // Add the sentinel.
  offsets_.push_back(current_offset);
  CHECK_GE(offsets_.back(), 0);
}

template <typename ColorType>
Ldi<ColorType>::Ldi(const ion::math::Vector2i& size,
                    std::vector<ColorType> colors, std::vector<float> depths)
    : size_(size), colors_(std::move(colors)), depths_(std::move(depths)) {
  CHECK_GE(size_[0], 0);
  CHECK_GE(size_[1], 0);
  CHECK_EQ(colors_.size(), size_[0] * size_[1]);
  CHECK_EQ(depths_.size(), size_[0] * size_[1]);
  offsets_.resize(size_[0] * size_[1] + 1);
  std::iota(offsets_.begin(), offsets_.end(), 0);
}

template <typename ColorType>
void Ldi<ColorType>::GetSampleCountAndOffset(const ion::math::Point2i& pixel,
                                             int* offset,
                                             int* sample_count) const {
  const int index = pixel[1] * size_[0] + pixel[0];
  *offset = offsets_[index];
  *sample_count = offsets_[index + 1] - *offset;
}

template <typename ColorType>
int Ldi<ColorType>::GetSampleCount(const ion::math::Point2i& pixel) const {
  int offset, sample_count;
  GetSampleCountAndOffset(pixel, &offset, &sample_count);
  return sample_count;
}

template <typename ColorType>
int Ldi<ColorType>::GetSampleCount() const {
  return offsets_.back();
}

template <typename ColorType>
absl::Span<const ColorType> Ldi<ColorType>::GetColors(
    const ion::math::Point2i& pixel) const {
  int offset, sample_count;
  GetSampleCountAndOffset(pixel, &offset, &sample_count);
  if (sample_count == 0) {
    return absl::Span<const ColorType>();
  } else {
    return absl::Span<const ColorType>(&colors_[offset], sample_count);
  }
}

template <typename ColorType>
absl::Span<ColorType> Ldi<ColorType>::GetMutableColors(
    const ion::math::Point2i& pixel) {
  int offset, sample_count;
  GetSampleCountAndOffset(pixel, &offset, &sample_count);
  if (sample_count == 0) {
    return absl::Span<ColorType>();
  } else {
    return absl::Span<ColorType>(&colors_[offset], sample_count);
  }
}

template <typename ColorType>
absl::Span<const float> Ldi<ColorType>::GetDepths(
    const ion::math::Point2i& pixel) const {
  int offset, sample_count;
  GetSampleCountAndOffset(pixel, &offset, &sample_count);
  if (sample_count == 0) {
    return absl::Span<const float>();
  } else {
    return absl::Span<const float>(&depths_[offset], sample_count);
  }
}

using Ldi1ui8 = Ldi<base::Color1ui8>;
using Ldi1f = Ldi<base::Color1f>;
using Ldi3ui8 = Ldi<base::Color3ui8>;
using Ldi3f = Ldi<base::Color3f>;
using Ldi4ui8 = Ldi<base::Color4ui8>;
using Ldi4f = Ldi<base::Color4f>;

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_LDI_H_
