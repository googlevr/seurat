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

#ifndef VR_SEURAT_GEOMETRY_KDTREE_H_
#define VR_SEURAT_GEOMETRY_KDTREE_H_

#include <memory>
#include <vector>

#include "ion/math/vector.h"
#include "absl/types/span.h"

namespace seurat {
namespace geometry {

// A wrapper around nanoflann for K-nearest neighbor and radius queries in a
// point cloud.
template <int Dimension>
class KdTree {
 public:
  using PointT = ion::math::Point<Dimension, float>;

  // Constructs a KdTree given the |points| array.
  explicit KdTree(absl::Span<const PointT> points);

  ~KdTree();

  // Returns in |result| the indices of the nearest neighbors to the
  // |query_point|. The size of the result is at most |query_size|. The
  // neighbors are _not_ sorted by distance to the query point.
  void KnnSearch(const PointT& query_point, int query_size,
                 std::vector<int>* result) const;

  // Returns in |result| the index of the nearest neighbor to the |query_point|.
  //
  // Returns true upon success, false if there are no points.
  //
  // Performance Note:  For single-neighbor queries, this is faster than
  //                    KnnSearch.
  bool NNSearch(const PointT& query_point, int* result) const;

  // Returns in |result| the indices of all points that have
  // a squared distance to the |query_point| which is less than
  // |search_radius_squared|. If |sorted| is true the points will be sorted by
  // distance to the query point.
  void RadiusSearch(const PointT& query_point, float search_radius_squared,
                    bool sorted, std::vector<int>* result) const;

  // Disallow copy and assign.
  KdTree(KdTree&&) = delete;
  KdTree& operator=(KdTree&&) = delete;

 private:
  // We use PIMPL here to screen KdTree clients from the compilation
  // dependencies of nanoflann, e.g. the requirement for exception handling.
  class Impl;

  std::unique_ptr<Impl> impl_;
};

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_KDTREE_H_
