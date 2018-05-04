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

#include "seurat/geometry/polygon.h"

#include <array>
#include <limits>

#include "ion/math/angle.h"
#include "ion/math/range.h"
#include "ion/math/rangeutils.h"
#include "ion/math/vector.h"
#include "seurat/base/color.h"
#include "seurat/base/util.h"
#include "seurat/geometry/convex_hull2d.h"
#include "seurat/image/image.h"

namespace seurat {
namespace geometry {

using ion::math::Anglef;
using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Range2f;
using ion::math::Vector2f;
using seurat::base::Color4f;
using seurat::image::Image4f;

Range2f Extent(absl::Span<const ion::math::Point2f> polygon,
               const std::array<Vector2f, 2>& direction, const Point2f& point) {
  Range2f extent;
  for (const auto& p : polygon) {
    Vector2f pos_vec = p - point;
    extent.ExtendByPoint(
        Point2f(Dot(pos_vec, direction[0]), Dot(pos_vec, direction[1])));
  }
  return extent;
}

float PolygonArea(absl::Span<const ion::math::Point2f> polygon) {
  // Use Stokes' Theorem to compute the area as $A = \cint x dy$.
  int vertex_count = polygon.size();
  float result = 0.0f;
  for (int curr = 0; curr < vertex_count; ++curr) {
    int next = (curr + 1) % vertex_count;
    result += 0.5f * (polygon[next][0] + polygon[curr][0]) *
              (polygon[next][1] - polygon[curr][1]);
  }
  return result;
}

Quad2f ComputeConvexPolygonOBB(absl::Span<const ion::math::Point2f> cpoly,
                               float* obb_area) {
  std::array<Vector2f, 2> best_dir;
  Point2f best_pnt;
  Range2f best_range;
  int cpoly_size = cpoly.size();
  CHECK_LE(3, cpoly_size);
  float best_box_area = std::numeric_limits<float>::infinity();
  for (int i = 0; i < cpoly_size; ++i) {
    const Point2f& pnt = cpoly[i];
    int next_i = (i + 1) % cpoly_size;
    // Compute extents in two orthogonal directions. The first direction is
    // along the current polygon edge.
    std::array<Vector2f, 2> dir;
    dir[0] = ion::math::Normalized(cpoly[next_i] - pnt);
    CHECK_NE(Vector2f::Zero(), dir[0]);
    // The second direction vector is obtained by rotating the first direction
    // vector counter-clockwise by 90 degrees.
    dir[1] = ion::math::Orthogonal(dir[0]);
    const Range2f& range = Extent(cpoly, dir, pnt);
    float box_area = ion::math::NVolume(range);
    if (box_area < best_box_area) {
      best_dir = dir;
      best_pnt = pnt;
      best_range = range;
      best_box_area = box_area;
    }
  }
  *obb_area = best_box_area;

  const Point2f& min = best_range.GetMinPoint();
  const Point2f& max = best_range.GetMaxPoint();

  Quad2f oriented_bbox;

  oriented_bbox[0] =
      Point2f(best_pnt + min[0] * best_dir[0] + min[1] * best_dir[1]);
  oriented_bbox[1] =
      Point2f(best_pnt + max[0] * best_dir[0] + min[1] * best_dir[1]);
  oriented_bbox[2] =
      Point2f(best_pnt + max[0] * best_dir[0] + max[1] * best_dir[1]);
  oriented_bbox[3] =
      Point2f(best_pnt + min[0] * best_dir[0] + max[1] * best_dir[1]);

  return oriented_bbox;
}

std::vector<Point2f> ComputeTextureBoundary(const Image4f& texture,
                                            int* opaque_texel_count) {
  // Build the list of points to be used as input to computing the 2d convex
  // hull. March from left to right adding the top and bottom point in every
  // non-empty column of texels. Construct the upper and lower chains of points.
  const int kLower = 0, kUpper = 1;
  std::array<std::vector<Point2f>, 2> chain;
  *opaque_texel_count = 0;
  for (int x = 0; x < texture.Width(); ++x) {
    int min_y = -1, max_y = -1;
    for (int y = 0; y < texture.Height(); ++y) {
      float alpha = texture.At(x, y)[3];
      if (alpha <= 0.0f) continue;
      ++*opaque_texel_count;
      if (min_y < 0) min_y = y;
      if (max_y < y) max_y = y;
    }
    // Add no points if all texels in this column are clear.
    if (min_y < 0) continue;

    // Add points at pixel's corners to ensure that the pixel will be fully
    // enclosed.

    chain[kLower].push_back(Point2f(x, min_y));
    chain[kLower].push_back(Point2f(x + 1, min_y));
    chain[kUpper].push_back(Point2f(x, max_y + 1));
    chain[kUpper].push_back(Point2f(x + 1, max_y + 1));
  }

  std::vector<Point2f> points(chain[kLower]);

  points.insert(points.end(), chain[kUpper].begin(), chain[kUpper].end());
  auto last = std::unique(points.begin(), points.end());
  points.erase(last, points.end());

  return points;
}

float ComputeOBBTilt(absl::Span<const Point2f> obb, int* pivot_index) {
  CHECK_EQ(4, obb.size());
  // Find either one of two longer edges of the oriented bounding box.
  Vector2f max_edge = obb[1] - obb[0];
  float max_len_sqr = LengthSquared(max_edge);
  const Vector2f edge = obb[2] - obb[1];
  float len_sqr = LengthSquared(edge);
  if (max_len_sqr < len_sqr) {
    max_edge = edge;
    max_len_sqr = len_sqr;
  }
  // Determine the angle phi in the interval [0, PI) so that a CCW rotation by
  // phi will make the oriented bounding box vertical with a portrait
  // orientation.
  float phi = 0.5f * M_PI - std::atan2(max_edge[1], max_edge[0]);
  if (phi < 0.0f) phi += M_PI;
  if (M_PI <= phi) phi -= M_PI;
  // Determine the index of the vertex of the oriented bounding box that will be
  // the bottom-left vertex after the rotation. For this, we rotate the bounding
  // box by the angle that we have just determined and see which of the vertices
  // of the bounding box ends up closest to the bottom-left corner of the
  // rotated bonding box.
  float cos_phi = std::cos(phi);
  float sin_phi = std::sin(phi);
  std::vector<Point2f> rot_obb;
  Range2f rot_range;
  for (const auto& p : obb) {
    Point2f rp(cos_phi * p[0] - sin_phi * p[1],
               sin_phi * p[0] + cos_phi * p[1]);
    rot_obb.push_back(rp);
    rot_range.ExtendByPoint(rp);
  }
  *pivot_index = -1;
  const Point2f& rot_min_pnt = rot_range.GetMinPoint();
  float min_dist_sqr = std::numeric_limits<float>::infinity();
  for (int i = 0; i < 4; ++i) {
    float dist_sqr = ion::math::LengthSquared(rot_obb[i] - rot_min_pnt);
    if (dist_sqr < min_dist_sqr) {
      *pivot_index = i;
      min_dist_sqr = dist_sqr;
    }
  }
  CHECK_GE(*pivot_index, 0);
  return phi;
}

bool IsConvexCounterClockwise(absl::Span<const ion::math::Point2f> points) {
  if (points.size() < 3) {
    return true;
  }
  for (int i = 0; i < points.size(); ++i) {
    const Point2f& a = points[i];
    const Point2f& b = points[(i + 1) % 4];
    const Point2f& c = points[(i + 2) % 4];
    Vector2f v1 = c - b;
    Vector2f v2 = a - b;
    // v1 'cross' v2 should have positive z-coordinate.
    if (v1[0] * v2[1] - v1[1] * v2[0] < 0.0f) {
      return false;
    }
  }
  return true;
}

}  // namespace geometry
}  // namespace seurat
