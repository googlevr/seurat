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

#include "seurat/baker/framework/ray_classifier.h"

#include "ion/math/transformutils.h"
#include "seurat/base/array2d.h"
#include "seurat/base/array2d_util.h"
#include "seurat/base/parallel.h"
#include "seurat/geometry/triangle.h"
#include "seurat/image/image.h"

namespace seurat {
namespace baker {

using base::Array2D;
using base::Color4f;
using base::MutableArray2DView;
using geometry::Plane3f;
using geometry::Quad3f;
using geometry::Raytracer;
using image::ImageView1f;
using image::ImageView4f;
using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Range1f;
using ion::math::Range2i;
using ion::math::Range3f;
using ion::math::Vector2i;
using ion::math::Vector3f;
using ion::math::Vector4f;

void ProjectingRayClassifier::Init(absl::Span<const Frame> frames) {
  frames_ = frames;

  // Build a mesh from all frames.
  std::vector<ion::math::Point3f> vertex_buffer;
  std::vector<int> index_buffer;
  for (const auto& frame : frames_) {
    const int base_index = vertex_buffer.size();
    vertex_buffer.insert(vertex_buffer.end(), frame.quad.begin(),
                         frame.quad.end());

    // Append two triangles:  012, 023
    //
    // 3----2
    // |   /|
    // |  / |
    // | /  |
    // 0----1

    index_buffer.push_back(base_index + 0);
    index_buffer.push_back(base_index + 1);
    index_buffer.push_back(base_index + 2);

    index_buffer.push_back(base_index + 0);
    index_buffer.push_back(base_index + 2);
    index_buffer.push_back(base_index + 3);
  }
  raytracer_ = Raytracer::Build(vertex_buffer, index_buffer);
}

namespace {

// Given a vector of samples-per-frame, merges all samples-per-frame vectors
// into a single output, with samples sorted lexicographically according to
// their pixel coordinates.
template <typename T>
void CombineSamplesPerFrame(
    int thread_count,
    const std::vector<std::vector<std::vector<T>>>& samples_per_frame_per_chunk,
    absl::Span<std::vector<T>> sorted_samples_per_frame) {
  const int num_frames = sorted_samples_per_frame.size();
  base::ParallelFor(thread_count, num_frames, [&](int frame_index) {
    std::vector<T>& all_samples = sorted_samples_per_frame[frame_index];
    for (const auto& partial_samples_per_frame : samples_per_frame_per_chunk) {
      const auto& partial_samples = partial_samples_per_frame[frame_index];
      all_samples.insert(all_samples.end(), partial_samples.begin(),
                         partial_samples.end());
    }
    std::sort(all_samples.begin(), all_samples.end());
  });
}

// Returns the normalized distance of a ray intersection from the ray's
// endpoint.
//
// In other words, it is the distance between the hit point & the endpoint,
// scaled such that a distance of 1 equals endpoint - origin units.
float NormalizedRayDistance(float t_hit) { return std::fabs(t_hit - 1.0f); }

}  // namespace

void ProjectingRayClassifier::CollectSolidSamples(
    const RayBundle& bundle,
    absl::Span<std::vector<RayBundle::RayIntersectionIndex>>
        solid_samples_per_frame,
    absl::Span<std::vector<int>> primary_frames_per_ray) const {
  const int num_rays = bundle.GetRayCount();
  const int num_frames = solid_samples_per_frame.size();

  // Parallelize by duplicating solid_samples_per_frame per-thread &
  // merge them afterwards.
  std::vector<std::vector<std::vector<RayBundle::RayIntersectionIndex>>>
      solid_samples_per_frame_per_thread(thread_count_);
  base::ParallelFor(thread_count_, thread_count_, [&](int tid) {
    std::vector<std::vector<RayBundle::RayIntersectionIndex>>&
        solid_samples_per_frame = solid_samples_per_frame_per_thread[tid];
    solid_samples_per_frame.resize(num_frames);
    std::vector<Raytracer::Intersection> intersections;
    for (int r = tid; r < num_rays; r += thread_count_) {
      const int num_intersections = bundle.GetIntersectionCount(r);
      primary_frames_per_ray[r].reserve(num_intersections);
      for (int i = 0; i < num_intersections; ++i) {
        Point3f origin = Point3f::Zero();
        Point3f endpoint = bundle.GetIntersectionPoint(r, i);

        raytracer_->FindAllIntersections(origin, endpoint - origin,
                                         &intersections);

        if (intersections.empty()) {
          continue;
        }

        // Assign to the primary frame.
        const Raytracer::Intersection primary =
            *std::min_element(intersections.begin(), intersections.end(),
                              [](const Raytracer::Intersection& lhs,
                                 const Raytracer::Intersection& rhs) {
                                return NormalizedRayDistance(lhs.t_hit) <
                                       NormalizedRayDistance(rhs.t_hit);
                              });
        const int primary_frame = primary.triangle_index / 2;
        primary_frames_per_ray[r].push_back(primary_frame);
        solid_samples_per_frame[primary_frame].push_back(std::make_tuple(r, i));

        // Assign to secondary frames.
        for (const auto& intersection : intersections) {
          int frame_index = intersection.triangle_index / 2;
          if (frame_index == primary_frame) {
            continue;
          }
          if (NormalizedRayDistance(intersection.t_hit) <
              secondary_frame_threshold_) {
            solid_samples_per_frame[frame_index].push_back(
                std::make_tuple(r, i));
          }
        }
      }
    }
  });

  CombineSamplesPerFrame(thread_count_, solid_samples_per_frame_per_thread,
                         solid_samples_per_frame);
}

void ProjectingRayClassifier::CollectFreespaceRays(
    const RayBundle& bundle,
    absl::Span<const std::vector<RayBundle::RayIntersectionIndex>>
        solid_rays_per_frame,
    absl::Span<const std::vector<int>> primary_frames_per_ray,
    absl::Span<std::vector<int>> freespace_rays_per_frame) const {
  const int num_rays = bundle.GetRayCount();
  const int num_frames = freespace_rays_per_frame.size();

  // Parallelize by duplicating freespace_rays_per_frame per-thread & merge them
  // afterwards.
  std::vector<std::vector<std::vector<int>>>
      freespace_rays_per_frame_per_thread(thread_count_);
  base::ParallelFor(thread_count_, thread_count_, [&](int tid) {
    std::vector<std::vector<int>>& freespace_rays_per_frame =
        freespace_rays_per_frame_per_thread[tid];
    freespace_rays_per_frame.resize(num_frames);
    std::vector<Raytracer::Intersection> intersections;
    for (int r = tid; r < num_rays; r += thread_count_) {
      Point3f origin = bundle.GetOrigin(r);
      if (bundle.GetIntersectionCount(r) == 0) {
        // There are no samples, so all intersections should be freespace
        // constraints.
        //
        // This is necessary to support processing partial scenes.  Pixels
        // without samples are assumed to be masked out and must carve
        // silhouettes through all geometry.
        raytracer_->FindAllIntersections(origin, bundle.GetDirection(r),
                                         &intersections);
        for (const auto& intersection : intersections) {
          int frame_index = intersection.triangle_index / 2;
          freespace_rays_per_frame[frame_index].push_back(r);
        }
        continue;
      }
      Point3f endpoint = bundle.GetIntersectionPoint(r, 0);
      raytracer_->FindAllIntersections(origin, endpoint - origin,
                                       &intersections);

      if (intersections.empty()) {
        continue;
      }

      for (const auto& intersection : intersections) {
        // TODO(puneetl):  Make this work for LDIs.
        int frame_index = intersection.triangle_index / 2;
        bool use_intersection = true;
        // Don't use this if it should be a "secondary sample", based on the
        // distance along the ray from the generating camera.
        //
        // Note that doing this check here, as opposed to only in
        // CollectSolidSamples(), allows for grazing-angle samples to remove
        // "fins" which may otherwise protrude from geometry.
        if (NormalizedRayDistance(intersection.t_hit) <
            secondary_frame_threshold_) {
          use_intersection = false;
        }
        if (std::binary_search(primary_frames_per_ray[r].begin(),
                               primary_frames_per_ray[r].end(), frame_index)) {
          // Skip this ray, since it is already assigned to this frame as a
          // primary solid sample.
          use_intersection = false;
        }
        if (!use_intersection) continue;
        for (int primary_frame_index : primary_frames_per_ray[r]) {
          CHECK_NE(frame_index, primary_frame_index);
          bool carve_silhouette = false;
          if (rendering_mode_ == RenderingMode::kZBuffer) {
            float t_hit_primary;
            // If the result will be rendered using a proper z-buffer, then
            // classify this ray as a freespace sample iff it hits before the
            // primary intersection.
            if (PlaneFromFrame(frames_[primary_frame_index])
                    .IntersectRay(origin, endpoint - origin, &t_hit_primary)) {
              carve_silhouette = intersection.t_hit < t_hit_primary;
            }
          } else {
            // If the result will *not* be rendered using a proper z-buffer,
            // then classify this ray as a freespace sample iff the quad renders
            // after the primary quad.
            carve_silhouette = frames_[frame_index].draw_order >
                               frames_[primary_frame_index].draw_order;
          }
          if (carve_silhouette) {
            freespace_rays_per_frame[frame_index].push_back(r);
            break;
          }
        }
      }
    }
  });

  CombineSamplesPerFrame(thread_count_, freespace_rays_per_frame_per_thread,
                         freespace_rays_per_frame);
}

std::vector<RayClassifier::ClassifiedRays>
ProjectingRayClassifier::ClassifyRays(const RayBundle& bundle) const {
  const int ray_count = bundle.GetRayCount();
  std::vector<std::vector<int>> primary_frames_per_ray(ray_count);

  std::vector<std::vector<RayBundle::RayIntersectionIndex>>
      solid_samples_per_frame(frames_.size());
  CollectSolidSamples(bundle, absl::MakeSpan(solid_samples_per_frame),
                      absl::MakeSpan(primary_frames_per_ray));

  std::vector<std::vector<int>> freespace_rays_per_frame(frames_.size());
  CollectFreespaceRays(bundle, solid_samples_per_frame, primary_frames_per_ray,
                       absl::MakeSpan(freespace_rays_per_frame));

  std::vector<ClassifiedRays> classified_rays_per_frame(frames_.size());
  base::BalancedParallelFor(
      thread_count_, frames_.size(), [&](int frame_index) {
        std::vector<RayBundle::RayIntersectionIndex>& solid_samples =
            solid_samples_per_frame[frame_index];
        std::vector<int>& freespace_rays =
            freespace_rays_per_frame[frame_index];

        // Remove all solid samples which were used as freespace rays.
        // TODO(puneetl):  This could be done more efficiently using a merge.
        auto solid_end = std::remove_if(
            solid_samples.begin(), solid_samples.end(),
            [&](const RayBundle::RayIntersectionIndex& solid) {
              int ray_index = std::get<0>(solid);
              return std::binary_search(freespace_rays.begin(),
                                        freespace_rays.end(), ray_index);
            });
        solid_samples.erase(solid_end, solid_samples.end());

        classified_rays_per_frame[frame_index] = {
            std::move(solid_samples),
            std::move(freespace_rays_per_frame[frame_index])};
      });
  return classified_rays_per_frame;
}

void DilatingRayClassifier::Init(absl::Span<const Frame> frames) {
  std::vector<Vector2i> sizes(frames.size());
  texture_sizer_->ComputeTextureSizes(frames, absl::MakeSpan(sizes));
  dilated_frames_.clear();
  dilated_frames_.insert(dilated_frames_.begin(), frames.begin(), frames.end());
  for (int i = 0; i < dilated_frames_.size(); ++i) {
    // Texture pixels lie on the vertices of the Frame, so this must compensate
    // by subtracting 1.
    //
    // Textures which are 1x1 must be special-cased.
    float scale_x = texture_filter_radius_ / std::max(1, (sizes[i][0] - 1));
    float scale_y = texture_filter_radius_ / std::max(1, (sizes[i][1] - 1));
    const Frame frame = dilated_frames_[i];
    Frame new_frame = frame;
    new_frame.quad[0] += (frame.quad[0] - frame.quad[1]) * scale_x +
                         (frame.quad[0] - frame.quad[3]) * scale_y;
    new_frame.quad[1] += (frame.quad[1] - frame.quad[0]) * scale_x +
                         (frame.quad[1] - frame.quad[2]) * scale_y;
    new_frame.quad[2] += (frame.quad[2] - frame.quad[3]) * scale_x +
                         (frame.quad[2] - frame.quad[1]) * scale_y;
    new_frame.quad[3] += (frame.quad[3] - frame.quad[2]) * scale_x +
                         (frame.quad[3] - frame.quad[0]) * scale_y;
    dilated_frames_[i] = new_frame;
  }
  delegate_->Init(dilated_frames_);
}

}  // namespace baker
}  // namespace seurat
