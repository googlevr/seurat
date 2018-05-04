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

#include "seurat/geometry/convex_hull2d.h"

#include <list>

#include "ion/math/vectorutils.h"

using ion::math::Point2f;
using ion::math::Vector2f;

namespace {

bool LeftTurn(const Point2f& a, const Point2f& b, const Point2f& c) {
  return ion::math::Cross(b - a, c - b) > 0.0f;
}

// Add a new point |p| to the tentative |hull|. After that, keep removing the
// next to last point, if needed, until the last 3 points form a left turn.
void AppendAndTrim(const Point2f& p, std::list<Point2f>* hull) {
  hull->push_back(p);
  bool left_turn = false;
  while (hull->size() > 2 && !left_turn) {
    // Form a triple with the last three points of the tentative hull and keep
    // deleting the middle point until the points make a right turn.
    auto ti = hull->end();
    const auto ci = --ti;
    const auto bi = --ti;
    const auto ai = --ti;
    left_turn = LeftTurn(*ai, *bi, *ci);
    if (!left_turn) {
      hull->erase(bi);
    }
  }
}

}  // namespace

namespace seurat {
namespace geometry {

std::vector<Point2f> ComputeConvexHull(const std::vector<Point2f>& points) {
  std::vector<Point2f> result;

  // The output is empty when the set of points is degenerate (e.g. collinear
  // points).
  if (points.size() < 3) return result;

  // Order the points (x,y)-lexicographically.

  struct {
    bool operator()(const Point2f& a, const Point2f& b) {
      if (a[0] < b[0]) return true;
      if (a[0] > b[0]) return false;
      return (a[1] < b[1]);
    }
  } lexicographic_less;

  std::vector<Point2f> sorted_points(points);
  std::sort(sorted_points.begin(), sorted_points.end(), lexicographic_less);

  // Several traversal orders are possible. Here we construct the convex hull by
  // moving counter-clockwise, leaving the interior of the convex hull on our
  // left side.

  // Construct the upper hull incrementally by walking from right to left.
  std::list<Point2f> upper_hull;
  for (auto pi = sorted_points.crbegin(); pi != sorted_points.crend(); ++pi) {
    AppendAndTrim(*pi, &upper_hull);
  }

  // Construct the lower hull incrementally by walking from left to right.
  std::list<Point2f> lower_hull;
  for (auto p : sorted_points) {
    AppendAndTrim(p, &lower_hull);
  }

  // Handle the degenerate case of collinear points. At this point in the
  // execution of the algorithm the endpoints of the upper and lower hulls are
  // duplicated: the  last point in the upper/lower hull is a copy of the first
  // point in the lower/upper hull, respectively. In the degenerate case both
  // the upper and lower hulls are reduced to one segment each, which means a
  // total of 4 points if we account for the duplicates.

  if (upper_hull.size() + lower_hull.size() < 5) return result;

  // Concatenate the upper and lower hulls to construct the result. Since the
  // last point in the upper/lower hull is the first point in the lower/upper
  // hull, respectively, do not output the last point in either list to avoid
  // duplicating points.
  for (auto pi = upper_hull.begin(); pi != --upper_hull.end(); ++pi)
    result.push_back(*pi);
  for (auto pi = lower_hull.begin(); pi != --lower_hull.end(); ++pi)
    result.push_back(*pi);

  return result;
}

}  // namespace geometry
}  // namespace seurat
