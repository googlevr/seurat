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

#include "seurat/baker/framework/texture_sizer.h"

#include "ion/math/vectorutils.h"
#include "seurat/base/util.h"
#include "seurat/geometry/quad.h"
#include "seurat/image/fixed_width_atlaser.h"

namespace seurat {
namespace baker {

using geometry::Quad3f;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2i;

void AreaTextureSizer::ComputeTextureSizes(
    absl::Span<const Frame> frames,
    absl::Span<ion::math::Vector2i> sizes) const {
  CHECK_EQ(frames.size(), sizes.size());
  for (int i = 0; i < frames.size(); ++i) {
    const Quad3f& quad = frames[i].quad;
    sizes[i][0] = std::max(ion::math::Length(quad[1] - quad[0]),
                           ion::math::Length(quad[3] - quad[2]));
    sizes[i][1] = std::max(ion::math::Length(quad[2] - quad[1]),
                           ion::math::Length(quad[3] - quad[0]));
  }
}

void ProjectedAreaTextureSizer::ComputeTextureSizes(
    absl::Span<const Frame> frames,
    absl::Span<ion::math::Vector2i> sizes) const {
  // The radius of the origin centered sphere which is scaled such that 1 unit
  // corresponds to 1 screen-space pixel.
  //
  // Compute this by circumference / (2 * pi).
  const float radius = (pixels_per_degree_ * 360.0f) / (2.0f * M_PI);
  // The area of the unit sphere in terms of square pixels at the target
  // resolution.
  const float square_pixels_per_sphere = 4.0f * M_PI * radius * radius;

  CHECK_EQ(frames.size(), sizes.size());
  for (int frame_index = 0; frame_index < frames.size(); ++frame_index) {
    const Quad3f& quad = frames[frame_index].quad;

    Quad3f quad_projected;
    for (int i = 0; i < 4; ++i) {
      quad_projected[i] =
          ion::math::Normalized(quad[i] - Point3f::Zero()) + Point3f::Zero();
    }

    // Estimate the area of the quad when projected onto the unit sphere,
    // normalized such that 1 is the area of the unit sphere.
    float estimated_area = 0.5f *
                           ion::math::Length(ion::math::Cross(
                               quad_projected[2] - quad_projected[0],
                               quad_projected[3] - quad_projected[1])) /
                           (4.0f * M_PI);
    // The target side-length.
    float side_length = std::sqrt(estimated_area * square_pixels_per_sphere);

    int side_length_rounded = static_cast<int>(side_length + 0.5f);

    sizes[frame_index] = Vector2i(side_length_rounded, side_length_rounded);
  }
}

BucketTextureSizer::BucketTextureSizer(std::unique_ptr<TextureSizer> delegate,
                                       const Vector2i& block_size)
    : delegate_(std::move(delegate)), block_size_(block_size) {}

void BucketTextureSizer::ComputeTextureSizes(
    absl::Span<const Frame> frames,
    absl::Span<ion::math::Vector2i> sizes) const {
  CHECK_EQ(frames.size(), sizes.size());
  delegate_->ComputeTextureSizes(frames, sizes);
  for (Vector2i& size : sizes) {
    // Handle an initial size of zero.
    size[0] = std::max(1, size[0]);
    size[1] = std::max(1, size[1]);
    // Bucket into the least multiple of the block size which is greater than
    // the initial texture size.
    size[0] =
        ((size[0] + block_size_[0] - 1) / block_size_[0]) * block_size_[0];
    size[1] =
        ((size[1] + block_size_[1] - 1) / block_size_[1]) * block_size_[1];
  }
}

void ConstrainedAtlasTextureSizer::ComputeTextureSizes(
    absl::Span<const Frame> frames,
    absl::Span<ion::math::Vector2i> sizes) const {
  CHECK_EQ(frames.size(), sizes.size());
  // Use 100 iterations of bisection.
  const int kBisectionIterations = 100;

  std::vector<Point2i> layout(frames.size());

  // Bisect to find the best lower bound.
  float lower_bound = 0.0f;
  float upper_bound = 1.0f;
  for (int i = 0; i < kBisectionIterations; ++i) {
    double pivot = (lower_bound + upper_bound) / 2.0;

    multiscale_sizer_factory_(pivot)->ComputeTextureSizes(frames, sizes);
    Vector2i atlas_size_target = atlaser_->GetAtlasSizeTarget();
    Vector2i total_texture_size;
    atlaser_->LayoutTiles(sizes, &total_texture_size, absl::MakeSpan(layout));
    if ((total_texture_size[0] <= atlas_size_target[0] &&
         total_texture_size[1] <= atlas_size_target[1])) {
      lower_bound = pivot;
    } else {
      upper_bound = pivot;
    }
  }

  multiscale_sizer_factory_(lower_bound)->ComputeTextureSizes(frames, sizes);
}

}  // namespace baker
}  // namespace seurat
