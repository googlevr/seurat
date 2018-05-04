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

#ifndef VR_SEURAT_TILER_SUBDIVISION_H_
#define VR_SEURAT_TILER_SUBDIVISION_H_

#include <array>
#include <vector>

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/geometry/quad.h"
#include "seurat/tiler/point_set.h"

namespace seurat {
namespace tiler {

// Interface for a fixed tree-like recursive subdivision.
//
// There may be multiple root cells, and all nodes of the tree correspond to a
// subset of points.
class Subdivision {
 public:
  virtual ~Subdivision() = default;

  // Returns all root cells.
  virtual void GetRoots(std::vector<int>* roots) const = 0;

  // Returns all child cells of the specified |cell|.
  virtual void GetChildren(int cell, std::vector<int>* children) const = 0;

  // Returns the indices of points in the tree cell.
  virtual absl::Span<const int> GetPointsInCell(int cell) const = 0;

  // Returns a set of 4 corner ray vectors indicating the edges of the cell's
  // bounding frustum.
  //
  // Vectors are in counterclockwise order when looking away from the
  // origin.
  //
  // The returned vectors must be normalized, or all zero if the cell does not
  // have well-defined "rails" (e.g. for a root node which might contain all
  // points).
  virtual std::array<ion::math::Vector3f, 4> GetRails(int cell) const = 0;

  // Rebuilds the space partitioning, if necessary.
  //
  // GetPointsInCell is undefined before this is called.
  //
  // For the sake of caching allocations, the same instance may be reinitialized
  // and reused multiple times.
  virtual void Init(const PointSet& points) = 0;
};

// Organizes a PointSet into a pyramid by partitioning points based on their
// orientation relative to the origin.
//
// Points are projected onto the faces of an origin-centered cube-map, and then
// organized according to a uniform quadtree.
//
// So, there are 6 root nodes.  The depth of the subdivision is the depth of
// the quadtree in each cube face.
class CubemapQuadtreeSubdivision : public Subdivision {
 public:
  explicit CubemapQuadtreeSubdivision(int depth);

  // Subdivision implementation.
  void GetRoots(std::vector<int>* roots) const override;
  void GetChildren(int node, std::vector<int>* children) const override;
  absl::Span<const int> GetPointsInCell(int cell) const override;
  std::array<ion::math::Vector3f, 4> GetRails(int cell) const override;
  void Init(const PointSet& points) override;

 private:
  struct Node {
    // A Span into points_ containing all points.
    absl::Span<int> points;

    // The range within the cube face of points.
    //
    // A node for the root of a cube face has the range (-1,-1) to (1, 1).
    ion::math::Range2f range;

    // The cube face containing this node.
    int cube_face;

    // The index in nodes_ of the start of the sequence of child nodes.
    int children_begin;
    // The number of child nodes.
    int children_count;
  };

  // Partitions a set of point_indices, assumed to be within the same cube-face,
  // into the 4 quadrants relative to the specified pivot point.
  std::array<absl::Span<int>, 4> Partition(const PointSet& points,
                                           const ion::math::Point2f& pivot,
                                           int major_axis,
                                           absl::Span<int> point_indices) const;

  // Maximum depth of the tree.
  int depth_;

  // The id of the most-recent PointSet, used to invalidate cached data
  // structures.
  int point_set_id_;

  // The set of all point indices, partitioned according to the structure of the
  // quadtree.
  std::vector<int> points_;

  // All nodes of the subdivision.  The first 6 nodes are the roots of each cube
  // face.
  std::vector<Node> nodes_;
};

// Wraps another Subdivision to dilate the bounds of all cells by some amount.
class BoundsDilatingSubdivision : public Subdivision {
 public:
  explicit BoundsDilatingSubdivision(float dilation_radians,
                                     std::unique_ptr<Subdivision> delegate)
      : dilation_radians_(dilation_radians), delegate_(std::move(delegate)) {}

  // Subdivision implementation.
  void GetRoots(std::vector<int>* roots) const override {
    delegate_->GetRoots(roots);
  }
  void GetChildren(int node, std::vector<int>* children) const override {
    delegate_->GetChildren(node, children);
  }
  absl::Span<const int> GetPointsInCell(int cell) const override {
    return delegate_->GetPointsInCell(cell);
  }
  std::array<ion::math::Vector3f, 4> GetRails(int cell) const override;
  void Init(const PointSet& points) override { delegate_->Init(points); }

 private:
  // The amount to dilate cell bounds, in radians.
  const float dilation_radians_;

  // The Subdivision to dilate.
  const std::unique_ptr<Subdivision> delegate_;
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SUBDIVISION_H_
