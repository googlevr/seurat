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

#include "seurat/baker/framework/frame.h"

#include <algorithm>
#include <numeric>

#include "ion/gfx/image.h"
#include "ion/math/range.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/base/util.h"
#include "seurat/geometry/convex_hull2d.h"
#include "seurat/geometry/cube_face.h"
#include "seurat/geometry/plane.h"
#include "seurat/geometry/polygon.h"
#include "seurat/geometry/triangle.h"

namespace seurat {
namespace baker {

using geometry::Plane2f;
using geometry::Plane3f;
using geometry::Quad2f;
using geometry::Quad3f;
using geometry::Triangle2f;
using geometry::Triangle3f;
using ion::math::Matrix2f;
using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Range1f;
using ion::math::Range2f;
using ion::math::Vector2f;
using ion::math::Vector2i;
using ion::math::Vector3f;

namespace {

Matrix4f WorldFromTangent(const Plane3f& plane, Point3f point_on_plane) {
  const Vector3f normal = plane.GetNormal();
  const Vector3f tangent = plane.GetTangent();
  const Vector3f cotangent =
      ion::math::Normalized(ion::math::Cross(normal, tangent));

  Matrix4f world_from_tangent(
      tangent[0], cotangent[0], normal[0], point_on_plane[0],  //
      tangent[1], cotangent[1], normal[1], point_on_plane[1],  //
      tangent[2], cotangent[2], normal[2], point_on_plane[2],  //
      0.0f, 0.0f, 0.0f, 1.0f);                                 //
  return world_from_tangent;
}

}  // namespace

bool operator==(const Frame& a, const Frame& b) {
  return a.draw_order == b.draw_order && a.quad == b.quad;
}

void InitializeApproximateDrawOrder(absl::Span<Frame> frames) {
  // Frames closer to the origin have higher draw order (i.e. are drawn later).
  struct FrameCenter {
    int frame_index;
    float squared_distance_to_origin;
  };
  std::vector<FrameCenter> frame_centers(frames.size());
  for (int frame_index = 0; frame_index < frames.size(); ++frame_index) {
    const Frame& frame = frames[frame_index];
    Point3f mean_point =
        std::accumulate(frame.quad.begin(), frame.quad.end(), Point3f::Zero()) /
        4.0f;
    float squared_distance_to_origin =
        ion::math::LengthSquared(mean_point - Point3f::Zero());
    frame_centers[frame_index] = {frame_index, squared_distance_to_origin};
  }
  std::sort(frame_centers.begin(), frame_centers.end(),
            [](const FrameCenter& lhs, const FrameCenter& rhs) {
              // Note that using '>' results in sorting in *descending order*.
              return lhs.squared_distance_to_origin >
                     rhs.squared_distance_to_origin;
            });
  for (int center_index = 0; center_index < frame_centers.size();
       ++center_index) {
    frames[frame_centers[center_index].frame_index].draw_order = center_index;
  }
}

geometry::Plane3f PlaneFromFrame(const Frame& frame) {
  return geometry::PlaneFromTriangle<float>(
      {{frame.quad[0], frame.quad[1], frame.quad[2]}});
}

bool DilateFrame(float resolution, Frame* frame) {
  // The following admittedly uses a very loose heuristic to compute the
  // scale factor by which to transform the given frame.
  //
  // In practice, this works well enough since the precise amount by which our
  // proxy-geometry is dilated turns out to be relatively unimportant.
  //
  // Too little dilation results in noticeable cracks, but too much dilation
  // results in nearly-imperceptible changes in the final result.

  Point3f frame_center =
      std::accumulate(frame->quad.begin(), frame->quad.end(), Point3f::Zero()) /
      4.0f;
  Plane3f plane = PlaneFromFrame(*frame);

  const Point3f eye = Point3f::Zero();

  // Scale each vertex away from the quad center by an amount which should push
  // it approximately 1 screenspace pixel away.
  std::array<float, 4> scale_factors;
  for (int i = 0; i < 4; ++i) {
    float scale_factor =
        1.0f +
        ion::math::Length(frame->quad[i] - eye) /
            (ion::math::Length(frame->quad[i] - frame_center) * resolution);
    if (!std::isfinite(scale_factor)) return false;
    scale_factors[i] = scale_factor;
  }

  const Matrix4f world_from_tangent = WorldFromTangent(plane, frame_center);
  const Matrix4f tangent_from_world = ion::math::Inverse(world_from_tangent);

  Quad3f dilated_quad;
  for (int i = 0; i < 4; ++i) {
    float scale = scale_factors[i];
    Point3f p = ion::math::ProjectPoint(tangent_from_world, frame->quad[i]);
    p *= scale;
    dilated_quad[i] = ion::math::ProjectPoint(world_from_tangent, p);
  }
  frame->quad = dilated_quad;
  return true;
}

bool WorldToFrame(const Frame& frame, const Point3f& point_world,
                  Point3f* point_frame) {
  const std::array<float, 4> w = frame.texcoord_w;
  const Triangle3f tri1 = {{frame.quad[0], frame.quad[1], frame.quad[2]}};
  const Triangle3f tri1tex = {
      {{0.0f, 0.0f, w[0]}, {w[1], 0.0f, w[1]}, {w[2], w[2], w[2]}}};
  const Triangle3f tri2 = {{frame.quad[0], frame.quad[2], frame.quad[3]}};
  const Triangle3f tri2tex = {
      {{0.0f, 0.0f, w[0]}, {w[2], w[2], w[2]}, {0.0f, w[3], w[3]}}};

  const Point3f bary1 = geometry::BarycentricFromPoint(tri1, point_world);
  if (geometry::BarycentricCoordinatesInTriangle(bary1)) {
    *point_frame = geometry::PointFromBarycentric(tri1tex, bary1);
    return true;
  }

  const Point3f bary2 = geometry::BarycentricFromPoint(tri2, point_world);
  if (geometry::BarycentricCoordinatesInTriangle(bary2)) {
    *point_frame = geometry::PointFromBarycentric(tri2tex, bary2);
    return true;
  }

  return false;
}

bool FrameToWorld(const Frame& frame, const Point3f& point_frame,
                  Point3f* point_world) {
  const std::array<float, 4> w = frame.texcoord_w;
  const Triangle3f tri1 = {{frame.quad[0], frame.quad[1], frame.quad[2]}};
  const Triangle3f tri1tex = {
      {{0.0f, 0.0f, w[0]}, {w[1], 0.0f, w[1]}, {w[2], w[2], w[2]}}};
  const Triangle3f tri2 = {{frame.quad[0], frame.quad[2], frame.quad[3]}};
  const Triangle3f tri2tex = {
      {{0.0f, 0.0f, w[0]}, {w[2], w[2], w[2]}, {0.0f, w[3], w[3]}}};

  const Point3f bary1 = geometry::BarycentricFromPoint(tri1tex, point_frame);
  if (geometry::BarycentricCoordinatesInTriangle(bary1)) {
    *point_world = geometry::PointFromBarycentric(tri1, bary1);
    return true;
  }

  const Point3f bary2 = geometry::BarycentricFromPoint(tri2tex, point_frame);
  if (geometry::BarycentricCoordinatesInTriangle(bary2)) {
    *point_world = geometry::PointFromBarycentric(tri2, bary2);
    return true;
  }

  return false;
}

bool FreespaceRayToFrameSpace(const Frame& frame, const Point3f& start,
                              const Vector3f& direction, Point2f* frame_space) {
  float t_hit;
  const Plane3f plane = PlaneFromFrame(frame);
  if (!plane.IntersectRay(start, direction, &t_hit)) {
    return false;
  }
  const Point3f hit_point = direction * t_hit + start;
  Point3f frame_space_hom;
  if (!WorldToFrame(frame, hit_point, &frame_space_hom)) {
    return false;
  }
  *frame_space = Point2f(frame_space_hom[0] / frame_space_hom[2],
                         frame_space_hom[1] / frame_space_hom[2]);
  return true;
}

bool SolidRayToFrameSpace(const Frame& frame, const Point3f& start,
                          const Point3f& end, Point2f* frame_space) {
  // Note that |start| is not used in the current implementation, since we are
  // intersecting the ray from the origin to the ray's endpoint against the
  // Frame.
  const Point3f origin = Point3f::Zero();
  float t_hit;
  const Plane3f plane = PlaneFromFrame(frame);
  if (!plane.IntersectRay(origin, end - origin, &t_hit)) {
    return false;
  }
  const Point3f hit_point = (end - origin) * t_hit + origin;
  Point3f frame_space_hom;
  if (!WorldToFrame(frame, hit_point, &frame_space_hom)) {
    return false;
  }
  *frame_space = Point2f(frame_space_hom[0] / frame_space_hom[2],
                         frame_space_hom[1] / frame_space_hom[2]);
  return true;
}

}  // namespace baker
}  // namespace seurat
