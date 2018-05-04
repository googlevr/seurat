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

#include "seurat/baker/framework/frame_sorter.h"

#include <algorithm>
#include <numeric>

#include "ion/math/vector.h"
#include "seurat/base/array2d_util.h"
#include "seurat/base/parallel.h"
#include "seurat/geometry/raytracer.h"
#include "seurat/geometry/triangle.h"
#include "seurat/image/image.h"

namespace seurat {
namespace baker {

using geometry::Raytracer;
using ion::math::Point3f;
using ion::math::Vector3f;

namespace {

// Returns a raytracer which traces rays through the quads of the given
// |frames|.
//
// The i'th frame corresponds to the 2i and 2i+1 triangles.
std::unique_ptr<Raytracer> RaytracerForFrames(absl::Span<const Frame> frames) {
  std::vector<ion::math::Point3f> vertex_buffer;
  std::vector<int> index_buffer;
  for (const auto& frame : frames) {
    const int base_index = vertex_buffer.size();
    vertex_buffer.insert(vertex_buffer.end(), frame.quad.begin(),
                         frame.quad.end());
    index_buffer.push_back(base_index + 0);
    index_buffer.push_back(base_index + 1);
    index_buffer.push_back(base_index + 2);

    index_buffer.push_back(base_index + 0);
    index_buffer.push_back(base_index + 2);
    index_buffer.push_back(base_index + 3);
  }
  return Raytracer::Build(vertex_buffer, index_buffer);
}

// Given a point & direction vector, finds the *line* intersection which is
// closest to the point.
//
// Returns false if no hit was found.
//
// Note: t_hit may be negative to indicate an intersection which was behind the
// point.
bool FindClosestLineIntersection(const Raytracer& raytracer,
                                 const Point3f& point,
                                 const Vector3f& direction, float* t_hit,
                                 int* triangle_hit_index) {
  // The closest hit can be found by searching along the ray both "forward" and
  // "reverse".
  float t_max = std::numeric_limits<float>::infinity();
  float current_t_hit;
  int current_triangle_hit_index;
  bool success = false;
  if (raytracer.FindFirstHit(point, direction, t_max, &current_t_hit,
                             &current_triangle_hit_index)) {
    *t_hit = current_t_hit;
    *triangle_hit_index = current_triangle_hit_index;
    t_max = current_t_hit;
    success = true;
  }
  if (raytracer.FindFirstHit(point, -direction, t_max, &current_t_hit,
                             &current_triangle_hit_index)) {
    *t_hit = -current_t_hit;
    *triangle_hit_index = current_triangle_hit_index;
    success = true;
  }
  return success;
}

}  // namespace

void FrameSorter::ComputeDrawOrder(absl::Span<const Point3f> points,
                                   absl::Span<Frame> frames) {
  std::unique_ptr<Raytracer> raytracer = RaytracerForFrames(frames);

  // Build mappings between Frames and Points by tracing rays.
  std::vector<int> frame_per_point(points.size());
  base::BalancedParallelFor(thread_count_, points.size(), [&](int point_index) {
    const Point3f& p = points[point_index];
    float t_hit;
    int triangle_hit_index;
    if (!FindClosestLineIntersection(*raytracer, p, p - Point3f::Zero(), &t_hit,
                                     &triangle_hit_index) ||
        t_hit < -1.0f) {
      // If there is no intersection, or if the closest line intersection is on
      // the other side of the origin.
      frame_per_point[point_index] = -1;
      return;
    }
    int frame_index = triangle_hit_index / 2;
    frame_per_point[point_index] = frame_index;
  });
  std::vector<std::vector<int>> points_per_frame(frames.size());
  for (int point_index = 0; point_index < points.size(); ++point_index) {
    int frame_index = frame_per_point[point_index];
    if (frame_index == -1) continue;
    points_per_frame[frame_index].push_back(point_index);
  }

  // Compute the distance to the median point in each frame.
  struct FrameDistance {
    int frame_index;
    float distance;
    bool operator<(const FrameDistance& rhs) const {
      return distance < rhs.distance;
    }
  };
  std::vector<FrameDistance> distances(frames.size());
  base::BalancedParallelFor(thread_count_, frames.size(), [&](int frame_index) {
    std::vector<int>& point_indices = points_per_frame[frame_index];
    distances[frame_index].frame_index = frame_index;

    // If there are no points assigned to this frame, use the mean frame vertex
    // as a fallback.
    if (point_indices.empty()) {
      Point3f mean_point =
          std::accumulate(frames[frame_index].quad.begin(),
                          frames[frame_index].quad.end(), Point3f::Zero()) /
          4.0f;
      distances[frame_index].distance =
          ion::math::Length(mean_point - Point3f::Zero());
      return;
    }

    float sum_distances = 0.0f;
    for (int index : point_indices) {
      sum_distances += ion::math::Length(points[index] - Point3f::Zero());
    }
    float mean_distance = sum_distances / point_indices.size();
    distances[frame_index].distance = mean_distance;
  });
  std::sort(distances.begin(), distances.end());

  for (int i = 0; i < frames.size(); ++i) {
    // Draw the back frames first.
    int draw_order = (frames.size() - 1) - i;
    int frame_index = distances[i].frame_index;
    frames[frame_index].draw_order = draw_order;
  }
}

}  // namespace baker
}  // namespace seurat
