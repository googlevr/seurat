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

#include "seurat/tiler/subdivision.h"

#include <algorithm>
#include <numeric>

#include "ion/base/logging.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/base/parallel.h"
#include "seurat/base/util.h"
#include "seurat/geometry/plane.h"

namespace seurat {
namespace tiler {

using geometry::Plane3f;
using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Range1f;
using ion::math::Range1i;
using ion::math::Range2f;
using ion::math::Vector2f;
using ion::math::Vector3f;

CubemapQuadtreeSubdivision::CubemapQuadtreeSubdivision(int depth)
    : depth_(depth), point_set_id_(PointSet::kInvalidId) {
  CHECK_LT(depth, 15);
}

void CubemapQuadtreeSubdivision::GetRoots(std::vector<int>* roots) const {
  roots->clear();
  for (int i = 0; i < 6; ++i) {
    roots->push_back(i);
  }
}

void CubemapQuadtreeSubdivision::GetChildren(int node,
                                             std::vector<int>* children) const {
  const Node& n = nodes_[node];
  for (int i = 0; i < n.children_count; ++i) {
    children->push_back(n.children_begin + i);
  }
}

namespace {

int CubeFaceFromPoint(const Point3f& point) {
  int major_axis = base::MajorAxisFromPosition(point);
  return major_axis + 3 * (point[major_axis] >= 0.0f);
}

int MajorAxisFromCubeFace(int cube_face) { return cube_face % 3; }

int SignFromCubeFace(int cube_face) { return cube_face / 3 > 0 ? 1 : -1; }

}  // namespace

std::array<absl::Span<int>, 4> CubemapQuadtreeSubdivision::Partition(
    const PointSet& points, const ion::math::Point2f& pivot, int major_axis,
    absl::Span<int> point_indices) const {
  // First, partition according to 'x', resulting in:
  //  | begin()    middle_x    end() |
  //
  // Then, partition according to 'y', resulting in:
  //  | begin()    low_x_middle_y    middle_x    high_x_middle_y    end() |

  int x_axis = (major_axis + 1) % 3;
  int y_axis = (major_axis + 2) % 3;

  auto middle_x = std::partition(
      point_indices.begin(), point_indices.end(), [&](int point_index) {
        const Point3f& point = points.positions[point_index];
        return point[x_axis] / std::fabs(point[major_axis]) < pivot[0];
      });

  auto low_x_middle_y =
      std::partition(point_indices.begin(), middle_x, [&](int point_index) {
        const Point3f& point = points.positions[point_index];
        return point[y_axis] / std::fabs(point[major_axis]) < pivot[1];
      });

  auto high_x_middle_y =
      std::partition(middle_x, point_indices.end(), [&](int point_index) {
        const Point3f& point = points.positions[point_index];
        return point[y_axis] / std::fabs(point[major_axis]) < pivot[1];
      });

  return {{
      absl::Span<int>(
          high_x_middle_y,
          std::distance(high_x_middle_y, point_indices.end())),  // +x, +y
      absl::Span<int>(low_x_middle_y,
                      std::distance(low_x_middle_y, middle_x)),  // -x, +y
      absl::Span<int>(
          point_indices.begin(),
          std::distance(point_indices.begin(), low_x_middle_y)),  // -x, -y
      absl::Span<int>(middle_x,
                      std::distance(middle_x, high_x_middle_y))  // +x, -y
  }};
}

void CubemapQuadtreeSubdivision::Init(const PointSet& points) {
  if (points.id == point_set_id_) {
    return;
  }

  point_set_id_ = points.id;

  const int point_count = points.positions.size();

  nodes_.clear();
  points_.resize(point_count);
  std::iota(points_.begin(), points_.end(), 0);

  // First, partition according to cube-face.
  std::sort(points_.begin(), points_.end(), [&](int lhs, int rhs) {
    return CubeFaceFromPoint(points.positions[lhs]) <
           CubeFaceFromPoint(points.positions[rhs]);
  });

  // Ranges of indices into points_ containing the indices of the points in each
  // cube face.
  std::array<Range1i, 6> ranges_per_face;
  for (int i = 0; i < point_count; ++i) {
    int cube_face = CubeFaceFromPoint(points.positions[points_[i]]);
    ranges_per_face[cube_face].ExtendByPoint(i);
  }
  for (int face = 0; face < 6; ++face) {
    const Range1i& range = ranges_per_face[face];
    absl::Span<int> point_slice;
    if (!range.IsEmpty()) {
      point_slice =
          absl::Span<int>(&points_[range.GetMinPoint()], range.GetSize() + 1);
    }
    nodes_.push_back(
        {point_slice, Range2f({-1.0f, -1.0f}, {1.0f, 1.0f}), face, 0, 0});
  }

  struct SplitCandidate {
    // The node to split.
    int node;

    // The depth of this node (root node starts at 0).
    int depth;

    // A slice into points_ of the indices of the points in this node.
    absl::Span<int> points;

    // The center of the node to split.
    Point2f pivot;

    // The side-length of the cube-face range of the node to split.
    //
    // For example, the root node of a cube face has a size of 2, since it
    // spans from (-1, -1) to (1, 1).
    float size;
  };

  std::vector<SplitCandidate> worklist;
  // Add the 6 cube faces to the worklist to split.
  for (int i = 0; i < 6; ++i) {
    const Point2f cube_face_center = Point2f::Zero();
    worklist.push_back(
        SplitCandidate{i, 0, nodes_[i].points, cube_face_center, 2.0f});
  }

  while (!worklist.empty()) {
    SplitCandidate to_split = worklist.back();
    worklist.pop_back();
    if (to_split.depth >= depth_) continue;

    int cube_face = nodes_[to_split.node].cube_face;
    std::array<absl::Span<int>, 4> child_slices =
        Partition(points, to_split.pivot, MajorAxisFromCubeFace(cube_face),
                  to_split.points);

    float pivot_delta = to_split.size / 4.0f;
    std::array<Vector2f, 4> delta = {{
        Vector2f(pivot_delta, pivot_delta),    //
        Vector2f(-pivot_delta, pivot_delta),   //
        Vector2f(-pivot_delta, -pivot_delta),  //
        Vector2f(pivot_delta, -pivot_delta)    //
    }};
    const Range2f& parent_range = nodes_[to_split.node].range;
    const Vector2f range_size(to_split.size / 2.0f, to_split.size / 2.0f);
    const Point2f range_min = parent_range.GetMinPoint();
    const Point2f range_center = parent_range.GetCenter();
    std::array<Range2f, 4> ranges = {{
        Range2f::BuildWithSize(range_center, range_size),
        Range2f::BuildWithSize({range_min[0], range_center[1]}, range_size),
        Range2f::BuildWithSize(range_min, range_size),
        Range2f::BuildWithSize({range_center[0], range_min[1]},
                               range_size)  //
    }};
    int child_node_start_index = nodes_.size();
    for (int i = 0; i < 4; ++i) {
      nodes_.push_back({child_slices[i], ranges[i], cube_face, 0, 0});

      Point2f new_pivot = to_split.pivot + delta[i];
      worklist.push_back({child_node_start_index + i, to_split.depth + 1,
                          child_slices[i], new_pivot, to_split.size / 2.0f});
    }
    nodes_[to_split.node].children_begin = child_node_start_index;
    nodes_[to_split.node].children_count = 4;
  }
}

std::array<ion::math::Vector3f, 4> CubemapQuadtreeSubdivision::GetRails(
    int cell) const {
  Range2f range = nodes_[cell].range;
  int cube_face = nodes_[cell].cube_face;
  int major_axis = MajorAxisFromCubeFace(cube_face);
  int sign = SignFromCubeFace(cube_face);
  int x_axis = (major_axis + 1) % 3;
  int y_axis = (major_axis + 2) % 3;

  Vector3f direction;
  direction[major_axis] = sign;

  std::array<Vector3f, 4> rails;
  rails.fill(direction);
  rails[0][x_axis] = range.GetMinPoint()[0];
  rails[0][y_axis] = range.GetMinPoint()[1];
  rails[1][x_axis] = range.GetMaxPoint()[0];
  rails[1][y_axis] = range.GetMinPoint()[1];
  rails[2][x_axis] = range.GetMaxPoint()[0];
  rails[2][y_axis] = range.GetMaxPoint()[1];
  rails[3][x_axis] = range.GetMinPoint()[0];
  rails[3][y_axis] = range.GetMaxPoint()[1];
  for (Vector3f& dir : rails) {
    ion::math::Normalize(&dir);
  }

  Vector3f inside = ion::math::Cross(rails[1], rails[0]);
  Plane3f plane(Point3f::Zero(), inside);
  if (plane.SignedDistanceToPoint(rails[2] + Point3f::Zero()) < 0.0f) {
    std::reverse(rails.begin(), rails.end());
  }

  return rails;
}

absl::Span<const int> CubemapQuadtreeSubdivision::GetPointsInCell(
    int cell) const {
  return nodes_[cell].points;
}

std::array<ion::math::Vector3f, 4> BoundsDilatingSubdivision::GetRails(
    int cell) const {
  std::array<Vector3f, 4> original = delegate_->GetRails(cell);
  std::array<Vector3f, 4> dilated;
  for (int i = 0; i < 4; ++i) {
    Vector3f prev = original[(i - 1 + 4) % 4];
    Vector3f cur = original[i];
    Vector3f next = original[(i + 1) % 4];
    Vector3f prev_dir = ion::math::Normalized(prev - cur);
    Vector3f next_dir = ion::math::Normalized(next - cur);

    // This relies on a small-angle approximation.
    dilated[i] = ion::math::Normalized(cur - prev_dir * dilation_radians_ -
                                       next_dir * dilation_radians_);
  }
  return dilated;
}

}  // namespace tiler
}  // namespace seurat
