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

#include "seurat/tiler/build_partition.h"

#include "absl/types/span.h"
#include "seurat/base/parallel.h"

namespace seurat {
namespace tiler {

bool BuildPartition::operator==(const BuildPartition& rhs) const {
  BuildPartition canonicalized = *this;
  canonicalized.Canonicalize();

  BuildPartition canonicalized_rhs = rhs;
  canonicalized_rhs.Canonicalize();

  if (canonicalized.point_indices_ != canonicalized_rhs.point_indices_) {
    return false;
  }
  if (model_ != rhs.model_) {
    return false;
  }
  return true;
}

void BuildPartition::Canonicalize() {
  std::sort(point_indices_.begin(), point_indices_.end());
}

void BuildPartition::Clear() {
  const int kNoPoint = -1;

  point_indices_.clear();
  worst_fit_point_ = {kNoPoint, -std::numeric_limits<float>::infinity()};
  best_fit_point_ = {kNoPoint, std::numeric_limits<float>::infinity()};
  total_error_ = 0.0f;
}

void BuildPartition::AddPoint(int point_index, float error) {
  point_indices_.push_back(point_index);
  // To enforce deterministic results, regardless of the order points are
  // added, break ties by using point_index.
  if (error > worst_fit_point_.error ||
      (error == worst_fit_point_.error &&
       point_index > worst_fit_point_.point_index)) {
    worst_fit_point_ = {point_index, error};
  }
  if (error < best_fit_point_.error ||
      (error == best_fit_point_.error &&
       point_index < best_fit_point_.point_index)) {
    best_fit_point_ = {point_index, error};
  }
  total_error_ += error;
}

void CanonicalizePartitions(int thread_count,
                            absl::Span<BuildPartition> build_partitions) {
  base::ParallelFor(thread_count, build_partitions.size(), [&](int bp_index) {
    build_partitions[bp_index].Canonicalize();
  });
}

}  // namespace tiler
}  // namespace seurat
