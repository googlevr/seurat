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

#ifndef VR_SEURAT_TILER_BUILD_PARTITION_H_
#define VR_SEURAT_TILER_BUILD_PARTITION_H_

#include <utility>
#include <vector>

#include "absl/types/span.h"
#include "seurat/base/util.h"
#include "seurat/geometry/plane.h"
#include "seurat/tiler/geometry_model.h"

namespace seurat {
namespace tiler {

// An intermediate representation of a Partition storing extra metadata about
// the points assigned to it.
//
// Points are handled by their integer index in a PointSet.
class BuildPartition {
 public:
  BuildPartition() : model_() { Clear(); }

  explicit BuildPartition(const GeometryModel& model) : model_(model) {
    Clear();
  }

  // Comparison operators, mostly useful for testing purposes.
  //
  // Two BuildPartitions are equal if they have the same GeometryModel and sets
  // of points (order does not matter).
  bool operator!=(const BuildPartition& rhs) const { return !(*this == rhs); }
  bool operator==(const BuildPartition& rhs) const;

  // Puts the BuildPartition into canonical form by sorting all point-indices.
  //
  // This is useful for ensuring determinism if samples are added from different
  // threads.
  void Canonicalize();

  // Returns a pointer to the partition's geometry model.
  GeometryModel* GetMutableModel() { return &model_; }

  // Immutable access to the partition's geometry model.
  const GeometryModel& GetModel() const { return model_; }

  // Returns whether the partition is empty and has no points.
  bool Empty() const { return GetSize() == 0; }

  // Returns the number of points in the partition.
  int GetSize() const { return point_indices_.size(); }

  // Returns the index of the point with the worst (greatest) error, or -1.
  int GetWorstFitPoint() const { return worst_fit_point_.point_index; }

  // Returns the index of the point with the best (least) error, or -1.
  int GetBestFitPoint() const { return best_fit_point_.point_index; }

  // Returns a list of indices of points in this partition.
  const std::vector<int>& GetPointIndices() const { return point_indices_; }

  // Returns the total error of all points in this partition.
  float GetTotalError() const { return total_error_; }

  // Removes all points from this partition.
  void Clear();

  // Adds a point with the given |error| to this partition.
  void AddPoint(int point_index, float error);

 private:
  struct PointErrorPair {
    // The index of a point in a PointSet.
    int point_index;

    // The cost of assigning the point to this partition's GeometryModel.
    float error;
  };

  // The geometry model representing the points in this partition.
  GeometryModel model_;

  // Indices into a PointSet of the points in this partition.
  std::vector<int> point_indices_;

  // The (index, error) for the point with the highest error.
  //
  // The index is -1 if there are no points.
  PointErrorPair worst_fit_point_;

  // The (index, error) for the point with the lowest error.
  //
  // The index is -1 if there are no points.
  PointErrorPair best_fit_point_;

  // The total error for all points.
  float total_error_;
};

// Canonicalizes multiple partitions in parallel.
void CanonicalizePartitions(int thread_count,
                            absl::Span<BuildPartition> build_partitions);

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_BUILD_PARTITION_H_
