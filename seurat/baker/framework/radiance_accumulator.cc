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

#include "seurat/baker/framework/radiance_accumulator.h"

#include "seurat/base/array2d.h"
#include "seurat/image/inpainting.h"

namespace seurat {
namespace baker {

using base::Array2D;
using base::Color4d;
using base::Color4f;
using image::Image4f;
using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2f;
using ion::math::Vector2i;
using ion::math::Vector3f;

namespace {

// Returns the continuous texture space position of a point in frame space.
Point2f TextureFromFrame(const Point2f& point_frame,
                         const Vector2i& texture_size) {
  // A note on texture-coordinates;
  //
  // The frame is rasterized such that each vertex of the 3D quad has a
  // texel centered on it.
  //
  // Note that this representation cannot handle a 1x1 texture.  So, that is
  // special-cased to resample from the center of the frame.
  //
  // Continuous texture space ranges from [0.0, width) x [0.0, height). The
  // center of discrete pixel (x, y) is located at (x + 0.5, y + 0.5) in
  // continuous texture space.

  Point2f point_texture(0.5f, 0.5f);
  if (texture_size[0] > 1) {
    point_texture[0] = point_frame[0] * (texture_size[0] - 1) + 0.5f;
  }
  if (texture_size[1] > 1) {
    point_texture[1] = point_frame[1] * (texture_size[1] - 1) + 0.5f;
  }
  return point_texture;
}

}  // namespace

void RadianceAccumulator::Resolve(image::Image3f* rgb) const {
  rgb->Resize(rgbw_.GetSize());

  // Note that the alpha-channel is currently ignored.
  Image4f rgba(rgbw_.GetSize());
  Array2D<bool> mask(rgbw_.GetSize());
  for (int y = 0; y < rgba.Height(); ++y) {
    for (int x = 0; x < rgba.Width(); ++x) {
      double weight = rgbw_.At(x, y)[3];
      for (int c = 0; c < 3; ++c) {
        if (weight == 0.0) {
          // Mask this for inpainting.
          mask.At(x, y) = true;
          continue;
        }
        float value = rgbw_.At(x, y)[c] / weight;
        value = std::max(value, 0.0f);
        value = std::min(value, 1.0f);
        if (!std::isfinite(value)) {
          // Detect NaN and replace with the result of inpainting.
          mask.At(x, y) = true;
          continue;
        }
        rgba.At(x, y)[c] = value;
      }
    }
  }

  image::InpaintSmooth(mask, &rgba);

  // Copy out the result of inpainting.
  for (int y = 0; y < rgba.Height(); ++y) {
    for (int x = 0; x < rgba.Width(); ++x) {
      for (int c = 0; c < 3; ++c) {
        rgb->At(x, y)[c] = rgba.At(x, y)[c];
      }
    }
  }
}

void RadianceAccumulator::Add(
    const RayBundle& bundle,
    absl::Span<const RayBundle::RayIntersectionIndex> solid_samples,
    absl::Span<const int> freespace_rays) {
  for (const std::tuple<int, int>& ray_and_intersection : solid_samples) {
    const int ray_index = std::get<0>(ray_and_intersection);
    const int intersection_index = std::get<1>(ray_and_intersection);

    const Point3f ray_start = bundle.GetOrigin(ray_index);
    const Point3f ray_end =
        bundle.GetIntersectionPoint(ray_index, intersection_index);
    const Vector3f ray_dir = ray_end - ray_start;
    // The frame-space point of the ray's intersection.
    Point2f point_frame;
    if (!SolidRayToFrameSpace(*frame_, ray_start, ray_end, &point_frame)) {
      continue;
    }
    // The position of the ray's intersection in continuous texture space.
    Point2f point_texture = TextureFromFrame(point_frame, rgbw_.GetSize());

    // Distance from the sample's ray to the resampling eye.
    float radius_eye = ion::math::Length(
        (ray_start - eye_) + (ion::math::Dot(eye_ - ray_start, ray_dir) /
                              ion::math::LengthSquared(ray_dir) * ray_dir));

    // Compute the directional filter weight. This is used to preserve crisp
    // specular highlights, by giving higher weights to samples from cameras
    // that are close to |eye_|.
    //
    // If baking sharp specular highlights, this weight may become very small,
    // so use double-precision here.
    //
    // TODO(ernstm): This filter is dependent on scene scale. Evaluate the usage
    // of spherical Gaussians.
    double directional_weight =
        std::exp(static_cast<double>(-0.5 * (radius_eye / sigma2_eye_)));
    // Use a minimum weight to not output black if baking very sharp
    // highlights.
    const double kMinWeight = 1.0e-6;
    directional_weight = std::max(directional_weight, kMinWeight);

    // Compute the range of discrete pixels that contain the sample in their
    // filter support. The range is [x_min, x_max) x [y_min, y_max).
    //
    // Example: A point at x = 1.25, y = 1.25 and filter radius = 1.0 results in
    // x_min = 0, x_max = 2, y_min = 0, y_max = 2.
    //
    // Example: A point at x = 1.75, y = 1.75 and filter radius = 1.0 results in
    // x_min = 1, x_max = 3, y_min = 1, y_max = 3.
    const int x_min = std::round(point_texture[0] - pixel_filter_->GetRadius());
    const int x_max = std::round(point_texture[0] + pixel_filter_->GetRadius());
    const int y_min = std::round(point_texture[1] - pixel_filter_->GetRadius());
    const int y_max = std::round(point_texture[1] + pixel_filter_->GetRadius());

    // Fetch the sample's RGBA color.
    const Color4f& sample_color =
        bundle.GetIntersectionColor(ray_index, intersection_index);

    // Loop over the pixels that contain the sample in their filter support.
    for (int x = x_min; x < x_max; ++x) {
      for (int y = y_min; y < y_max; ++y) {
        const Point2i pixel(x, y);
        if (!rgbw_.IsInside(pixel)) {
          continue;
        }
        const Point2f pixel_center(x + 0.5f, y + 0.5f);
        const float pixel_weight = pixel_filter_->Eval(
            (point_texture - pixel_center) + Point2f::Zero());
        const double weight = directional_weight * pixel_weight;
        rgbw_.At(pixel) +=
            Color4d(sample_color[0], sample_color[1], sample_color[2], 1.0) *
            weight;
      }
    }
  }
}

}  // namespace baker
}  // namespace seurat
