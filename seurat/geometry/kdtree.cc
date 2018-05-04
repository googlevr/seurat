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

#include "seurat/geometry/kdtree.h"

#include "ion/math/vector.h"
#include "nanoflann.hpp"
#include "seurat/base/util.h"

namespace seurat {
namespace geometry {

namespace {

// An adaptor class to use Nanoflann with a vector of points.
template <int Dimension>
class PointCloudAdapter {
 public:
  using PointT = ion::math::Point<Dimension, float>;

  explicit PointCloudAdapter(absl::Span<const PointT> points)
      : points_(points) {}

  size_t kdtree_get_point_count() const { return points_.size(); }

  float kdtree_distance(const float* p1, const size_t idx_p2,
                        size_t size) const {
    float l2norm_sq = 0.0f;
    for (int d = 0; d < PointT::kDimension; ++d) {
      const float diff = p1[d] - points_.at(idx_p2)[d];
      l2norm_sq += diff * diff;
    }
    return l2norm_sq;
  }

  float kdtree_get_pt(const size_t idx, int dim) const {
    return points_.at(idx)[dim];
  }

  template <class BBOX>
  bool kdtree_get_bbox(const BBOX& bbox) const {
    // Returning false indicates that nanoflann should default to its standard
    // bbox method instead of any custom routine.
    return false;
  }

 private:
  const absl::Span<const PointT> points_;
};

}  // namespace

// Implementation class for KdTree. Defined here to screen clients of
// KdTree from the compilation dependencies of nanoflann, e.g. exception
// handling.
template <int Dimension>
class KdTree<Dimension>::Impl {
 public:
  // Constructs a KdTree for the cloud of |points|.
  explicit Impl(absl::Span<const PointT> points)
      : cloud_(points),
        index_adaptor_(PointT::kDimension, cloud_,
                       nanoflann::KDTreeSingleIndexAdaptorParams()) {
    index_adaptor_.buildIndex();
  }

  // Returns in |result| the indices of the nearest neighbors to the
  // |query_point|. The size of the result is at most |query_size|.
  void KnnSearch(const PointT& query_point, int query_size,
                 std::vector<int>* result) const;

  // Returns in |result| the index of the nearest neighbor to the
  // |query_point|.
  //
  // Returns true upon success, false if there is no query point.
  bool NNSearch(const PointT& query_point, int* result) const;

  // Returns in |result| the indices of all points that have
  // a squared distance to the |query_point| which is less than
  // |search_radius_squared|.
  void RadiusSearch(const PointT& query_point, float search_radius_squared,
                    bool sorted, std::vector<int>* result) const;

 private:
  // PointCloudAdapter object used by nanoflann as a data adaptor.
  PointCloudAdapter<Dimension> cloud_;

  // Nanoflann indexer object.
  nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<float, PointCloudAdapter<Dimension>>,
      PointCloudAdapter<Dimension>, PointT::kDimension>
      index_adaptor_;
};

// KdTree::Impl methods

template <int Dimension>
void KdTree<Dimension>::Impl::KnnSearch(const PointT& query_point,
                                        int query_size,
                                        std::vector<int>* result) const {
  const int point_count = cloud_.kdtree_get_point_count();
  // Nanoflann fails if initialized without any points, so we must handle that
  // edge-case separately.
  if (query_size >= point_count) {
    for (int i = 0; i < point_count; ++i) {
      result->push_back(i);
    }
    return;
  }
  const int sentinel = point_count;
  std::vector<size_t> neighbors(query_size, sentinel);
  std::vector<float> distances(query_size);
  index_adaptor_.knnSearch(query_point.Data(), query_size, neighbors.data(),
                           distances.data());
  for (const auto& neighbor : neighbors) {
    if (neighbor != sentinel) {
      result->push_back(neighbor);
    }
  }
}

template <int Dimension>
bool KdTree<Dimension>::Impl::NNSearch(const PointT& query_point,
                                       int* result) const {
  const int point_count = cloud_.kdtree_get_point_count();
  // Nanoflann fails if initialized without any points, so we must handle that
  // edge-case separately.
  if (point_count == 0) {
    return false;
  }
  size_t neighbor;
  float distance;
  index_adaptor_.knnSearch(query_point.Data(), 1, &neighbor, &distance);
  *result = neighbor;
  return true;
}

template <int Dimension>
void KdTree<Dimension>::Impl::RadiusSearch(const PointT& query_point,
                                           float search_radius_squared,
                                           bool sorted,
                                           std::vector<int>* result) const {
  std::vector<std::pair<size_t, float>> index_distance_pairs;
  nanoflann::SearchParams search_parameters;
  search_parameters.sorted = sorted;
  index_adaptor_.radiusSearch(query_point.Data(), search_radius_squared,
                              index_distance_pairs, search_parameters);
  for (const auto& item : index_distance_pairs) {
    result->push_back(item.first);
  }
}

// KdTree methods

template <int Dimension>
KdTree<Dimension>::KdTree(absl::Span<const PointT> points)
    : impl_(new KdTree::Impl(points)) {}

template <int Dimension>
KdTree<Dimension>::~KdTree() = default;

template <int Dimension>
void KdTree<Dimension>::KnnSearch(const PointT& query_point, int query_size,
                                  std::vector<int>* result) const {
  impl_->KnnSearch(query_point, query_size, result);
}

template <int Dimension>
bool KdTree<Dimension>::NNSearch(const PointT& query_point, int* result) const {
  return impl_->NNSearch(query_point, result);
}

template <int Dimension>
void KdTree<Dimension>::RadiusSearch(const PointT& query_point,
                                     float search_radius_squared, bool sorted,
                                     std::vector<int>* result) const {
  impl_->RadiusSearch(query_point, search_radius_squared, sorted, result);
}

// Explicit template instantiation.
template class KdTree<2>;
template class KdTree<3>;

}  // namespace geometry
}  // namespace seurat
