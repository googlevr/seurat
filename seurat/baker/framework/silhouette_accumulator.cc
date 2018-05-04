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

#include "seurat/baker/framework/silhouette_accumulator.h"

namespace seurat {
namespace baker {

using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Vector2f;
using ion::math::Vector2i;

void SilhouetteAccumulator::Resolve(image::Image1f* alpha) const {
  const std::unique_ptr<const ImplicitSilhouette> silhouette(
      silhouette_buffer_->Resolve().release());
  const int sample_count = supersample_factor_ * supersample_factor_;
  for (int y = 0; y < alpha->Height(); ++y) {
    for (int x = 0; x < alpha->Width(); ++x) {
      // Sample the implicit silhouette, counting the number of samples which
      // are solid.
      int solid_sample_count = 0;
      for (int v = 0; v < supersample_factor_; ++v) {
        for (int u = 0; u < supersample_factor_; ++u) {
          Point2f point = Point2f::Zero() +
                          Vector2f(x * supersample_factor_ + u + 0.5f,
                                   y * supersample_factor_ + v + 0.5f) /
                              Vector2f(alpha->GetSize() * supersample_factor_);
          if (silhouette->IsSolidAtPoint(point)) {
            solid_sample_count++;
          }
        }
      }
      alpha->At(x, y)[0] =
          static_cast<float>(solid_sample_count) / sample_count;
    }
  }
}

void SilhouetteAccumulator::Add(
    const RayBundle& bundle,
    absl::Span<const std::tuple<int, int>> solid_samples,
    absl::Span<const int> freespace_rays) {
  // Convert rays from the |bundle| into 2D frame-space samples & add them to
  // the buffer.

  // Handle solid rays.
  for (const auto& solid : solid_samples) {
    int ray_index = std::get<0>(solid);
    int intersection_index = std::get<1>(solid);

    const Point3f ray_start = bundle.GetOrigin(ray_index);
    const Point3f ray_end =
        bundle.GetIntersectionPoint(ray_index, intersection_index);
    Point2f point;
    if (!SolidRayToFrameSpace(*frame_, ray_start, ray_end, &point)) {
      continue;
    }
    silhouette_buffer_->AddSolidSample(point);
  }

  // Handle freespace rays.
  for (const int ray_index : freespace_rays) {
    const Point3f ray_start = bundle.GetOrigin(ray_index);

    Point2f point;
    if (!FreespaceRayToFrameSpace(*frame_, ray_start,
                                  bundle.GetDirection(ray_index), &point)) {
      continue;
    }
    silhouette_buffer_->AddFreespaceSample(point);
  }
}

}  // namespace baker
}  // namespace seurat
