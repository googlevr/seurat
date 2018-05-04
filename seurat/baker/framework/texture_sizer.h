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

#ifndef VR_SEURAT_BAKER_FRAMEWORK_TEXTURE_SIZER_H_
#define VR_SEURAT_BAKER_FRAMEWORK_TEXTURE_SIZER_H_

#include <functional>
#include <memory>
#include <vector>

#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/baker/framework/frame.h"
#include "seurat/image/atlaser.h"

namespace seurat {
namespace baker {

// Given a set of Frames, determines how much memory to allocate to their
// textures.
//
// Note that the texture size is computed "blind", with no knowledge of the
// texture's contents.  Later pipeline stages may perform further downsampling
// based on actual texture data.
class TextureSizer {
 public:
  virtual ~TextureSizer() = default;

  // Computes the size of the textures to use for each of the given |frames|.
  virtual void ComputeTextureSizes(
      absl::Span<const Frame> frames,
      absl::Span<ion::math::Vector2i> sizes) const = 0;
};

// Uses the area of each Frame's quad to estimate the texture size.
//
// Note:  This is not very practical for use on real-world data, but is useful
// for testing.
class AreaTextureSizer : public TextureSizer {
 public:
  AreaTextureSizer() = default;
  ~AreaTextureSizer() override = default;

  // TextureSizer implementation.
  void ComputeTextureSizes(
      absl::Span<const Frame> frames,
      absl::Span<ion::math::Vector2i> sizes) const override;
};

// Uses the size of each Frame, when projected onto the unit-sphere, to
// determine their texture-size.
class ProjectedAreaTextureSizer : public TextureSizer {
 public:
  explicit ProjectedAreaTextureSizer(float pixels_per_degree)
      : pixels_per_degree_(pixels_per_degree) {}
  ~ProjectedAreaTextureSizer() override = default;

  // TextureSizer implementation.
  void ComputeTextureSizes(
      absl::Span<const Frame> frames,
      absl::Span<ion::math::Vector2i> sizes) const override;

 private:
  const float pixels_per_degree_;
};

// Rounds up from the result of another TextureSizer into the nearest bucket
// size.
class BucketTextureSizer : public TextureSizer {
 public:
  BucketTextureSizer(std::unique_ptr<TextureSizer> delegate,
                     const ion::math::Vector2i& block_size);

  // TextureSizer implementation.
  void ComputeTextureSizes(
      absl::Span<const Frame> frames,
      absl::Span<ion::math::Vector2i> sizes) const override;

 private:
  // Used to compute initial texture sizes to be rounded.
  std::unique_ptr<TextureSizer> delegate_;

  // All bucket sizes are multiples of this block size.
  const ion::math::Vector2i block_size_;
};

// Bisects over possible TextureSizers until one is found which atlases to a
// total texture within the target size.
class ConstrainedAtlasTextureSizer : public TextureSizer {
 public:
  // Constructs a TextureSizer which must scale its output according to the
  // given |scale|.
  //
  // No precise "scaling" behavior is required, except that the outputs of the
  // texture sizer should be monotonically-increasing with respect to the
  // specified |scale|.
  using MultiscaleSizerFactory =
      std::function<std::unique_ptr<TextureSizer>(float scale)>;

  ConstrainedAtlasTextureSizer(std::shared_ptr<image::Atlaser> atlaser,
                               MultiscaleSizerFactory multiscale_sizer_factory)
      : atlaser_(std::move(atlaser)),
        multiscale_sizer_factory_(std::move(multiscale_sizer_factory)) {}

  // TextureSizer implementation.
  void ComputeTextureSizes(
      absl::Span<const Frame> frames,
      absl::Span<ion::math::Vector2i> sizes) const override;

 private:
  // Defines the process of laying out texture tiles in a texture atlas,
  // including possible constraints on the atlas size.
  std::shared_ptr<image::Atlaser> atlaser_;

  // Constructs TextureSizers parameterized by a scale factor.
  MultiscaleSizerFactory multiscale_sizer_factory_;
};

}  // namespace baker
}  // namespace seurat

#endif  // VR_SEURAT_BAKER_FRAMEWORK_TEXTURE_SIZER_H_
