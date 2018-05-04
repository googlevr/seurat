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

#include "seurat/tiler/partitioner_stages.h"

#include <algorithm>
#include <mutex>  // NOLINT(build/c++11)
#include <numeric>
#include <queue>

#include "seurat/base/parallel.h"
#include "seurat/base/reporting.h"
#include "seurat/base/util.h"
#include "seurat/geometry/kdtree.h"

namespace seurat {
namespace tiler {

using ion::math::Point3f;
using ion::math::Vector3f;

void RandomizedPartitionInitializationStage::Run(
    const PointSet& points, absl::Span<BuildPartition> build_partitions) {
  // Randomly select samples to form initial "seed" partitions.
  const int num_points = points.positions.size();
  std::vector<int> all_point_indices(num_points);
  std::iota(all_point_indices.begin(), all_point_indices.end(), 0);
  std::shuffle(all_point_indices.begin(), all_point_indices.end(), random_);

  const int partition_count =
      std::min(static_cast<int>(build_partitions.size()), num_points);

  base::ParallelFor(thread_count_, partition_count, [&](int bp_index) {
    const int point_index = all_point_indices[bp_index];
    BuildPartition& bp = build_partitions[bp_index];
    geometry_solver_->InitializeModel(point_index, bp.GetMutableModel());
  });
}

void GeometryModelRefinementStage::Run(
    const PointSet& points, absl::Span<BuildPartition> build_partitions) {
  // Parallelize over all partitions.
  base::ParallelFor(thread_count_, build_partitions.size(), [&](int bp_index) {
    BuildPartition& p = build_partitions[bp_index];
    if (p.Empty()) {
      return;
    }

    const std::vector<int>& points_in_partition = p.GetPointIndices();
    GeometryModel* model = p.GetMutableModel();
    bool success = geometry_solver_->FitModel(points_in_partition, model);
    if (!success) {
      geometry_solver_->InitializeModel(p.GetPointIndices().front(), model);
    }

    std::vector<int> point_indices = p.GetPointIndices();
    p.Clear();

    // Re-add the samples to the partition with updated error.
    for (const int point_index : point_indices) {
      float err = geometry_solver_->ComputeError(point_index, p.GetModel());
      p.AddPoint(point_index, err);
    }
  });

  CanonicalizePartitions(thread_count_, build_partitions);
}

void GreedyPointAssignmentStage::Run(
    const PointSet& point_set, absl::Span<BuildPartition> build_partitions) {
  std::vector<Point3f> normalized_partition_centers;
  normalized_partition_centers.reserve(build_partitions.size());
  for (const auto& p : build_partitions) {
    normalized_partition_centers.push_back(
        ion::math::Normalized(p.GetModel().center - Point3f::Zero()) +
        Point3f::Zero());
  }

  geometry::KdTree<3> kdtree(normalized_partition_centers);

  for (auto& bp : build_partitions) {
    bp.Clear();
  }

  std::vector<std::mutex> build_partition_mutex(build_partitions.size());
  base::ParallelFor(thread_count_, thread_count_, [&](int tid) {
    std::vector<int> partition_indices;
    for (int point_index = tid; point_index < point_set.positions.size();
         point_index += thread_count_) {
      partition_indices.clear();
      Point3f query_point =
          ion::math::Normalized(point_set.positions[point_index] -
                                Point3f::Zero()) +
          Point3f::Zero();
      kdtree.KnnSearch(query_point, neighboring_partition_count_,
                       &partition_indices);

      float best_error = std::numeric_limits<float>::max();
      const int kNoBestPartitionIdx = -1;
      int best_partition_idx = kNoBestPartitionIdx;
      for (int partition_index : partition_indices) {
        const GeometryModel& model =
            build_partitions[partition_index].GetModel();
        const float error = geometry_solver_->ComputeError(point_index, model);
        if (error < best_error) {
          best_error = error;
          best_partition_idx = partition_index;
        }
      }

      if (best_partition_idx == kNoBestPartitionIdx) {
        // Exit the lambda without assigning a partition on this iteration.
        // The point will have another chance if the partitioner has time to
        // run another iteration.
        //
        // TODO(puneetl) - b/28677934 - Decide what to do if the last
        // iteration doesn't assign a partition.
        return;
      }
      build_partition_mutex[best_partition_idx].lock();
      build_partitions[best_partition_idx].AddPoint(point_index, best_error);
      build_partition_mutex[best_partition_idx].unlock();
    }
  });

  CanonicalizePartitions(thread_count_, build_partitions);
}

namespace {

// Splits the given partition, modifying the |original| partition and
// returning the other new partition.
//
// Returns true if the split was successful, false otherwise.
bool Split(const GeometrySolver& geometry_solver, BuildPartition* original,
           BuildPartition* new_partition) {
  CHECK_NOTNULL(original);
  CHECK_NE(original->GetBestFitPoint(), original->GetWorstFitPoint());

  // Use the point with highest error to seed the new partition.
  GeometryModel* new_model = new_partition->GetMutableModel();
  *new_model = original->GetModel();
  geometry_solver.InitializeModel(original->GetWorstFitPoint(), new_model);
  new_partition->Clear();

  // The samples in the |original| partition to split.
  std::vector<int> original_points = original->GetPointIndices();

  original->Clear();

  // Assign each point to the better partition, according to the cost function.
  for (int point_index : original_points) {
    float error_original =
        geometry_solver.ComputeError(point_index, original->GetModel());
    float error_new = geometry_solver.ComputeError(point_index, *new_model);
    if (error_new < error_original) {
      new_partition->AddPoint(point_index, error_new);
    } else {
      original->AddPoint(point_index, error_original);
    }
  }

  return true;
}

}  // namespace

void PartitionSplittingStage::Run(const PointSet& point_set,
                                  absl::Span<BuildPartition> build_partitions) {
  struct CandidatePartitionToSplit {
    float total_error;
    BuildPartition* build_partition;
    bool operator<(const CandidatePartitionToSplit& rhs) const {
      return total_error < rhs.total_error;
    }
  };
  std::priority_queue<CandidatePartitionToSplit> candidate_heap;
  std::vector<BuildPartition*> dead_partitions;

  for (BuildPartition& p : build_partitions) {
    if (p.Empty()) {
      dead_partitions.push_back(&p);
    } else {
      candidate_heap.push({p.GetTotalError(), &p});
    }
  }

  int new_partitions = 0;
  for (BuildPartition* to_revive : dead_partitions) {
    if (candidate_heap.empty()) {
      return;
    }
    CandidatePartitionToSplit candidate_to_split = candidate_heap.top();
    candidate_heap.pop();

    if (candidate_to_split.build_partition->GetSize() <
        minimum_points_per_partition_) {
      // Too small to split.
      continue;
    }

    bool success =
        Split(*geometry_solver_, candidate_to_split.build_partition, to_revive);
    if (success) {
      candidate_heap.push({to_revive->GetTotalError(), to_revive});
      candidate_heap.push({candidate_to_split.build_partition->GetTotalError(),
                           candidate_to_split.build_partition});
      new_partitions++;
    }
  }
}

int PointExchangeStage::FindBestPartition(
    absl::Span<const BuildPartition> partitions, int point_index,
    float* cost) const {
  float best_error = std::numeric_limits<float>::infinity();
  int best_partition_index = 0;
  for (int partition_index = 0; partition_index < partitions.size();
       ++partition_index) {
    const GeometryModel& model = partitions[partition_index].GetModel();
    const float error = geometry_solver_->ComputeError(point_index, model);
    if (error <= best_error) {
      best_error = error;
      best_partition_index = partition_index;
    }
  }
  *cost = best_error;
  return best_partition_index;
}

void PointExchangeStage::Run(const PointSet& point_set,
                             absl::Span<BuildPartition> build_partitions) {
  relevant_point_indices_.clear();

  for (auto& p : build_partitions) {
    std::copy(p.GetPointIndices().begin(), p.GetPointIndices().end(),
              std::back_inserter(relevant_point_indices_));

    p.Clear();
  }

  std::vector<std::mutex> build_partition_mutex(build_partitions.size());
  base::ParallelFor(thread_count_, relevant_point_indices_.size(),
                    [&](int point_index_index) {
                      int point_index =
                          relevant_point_indices_[point_index_index];

                      float best_partition_cost;
                      int best_partition_index = FindBestPartition(
                          build_partitions, point_index, &best_partition_cost);

                      build_partition_mutex[best_partition_index].lock();
                      build_partitions[best_partition_index].AddPoint(
                          point_index, best_partition_cost);
                      build_partition_mutex[best_partition_index].unlock();
                    });

  CanonicalizePartitions(thread_count_, build_partitions);
}

void DepthBasedRedistributionStage::Run(
    const PointSet& point_set, absl::Span<BuildPartition> build_partitions) {
  relevant_point_indices_.clear();

  for (auto& p : build_partitions) {
    std::copy(p.GetPointIndices().begin(), p.GetPointIndices().end(),
              std::back_inserter(relevant_point_indices_));

    p.Clear();
  }

  const int point_count = relevant_point_indices_.size();
  const int partition_count = build_partitions.size();

  if (point_count == 0 || partition_count == 0) {
    return;
  }

  // Sort point indices based on depth to later assign points to partitions.
  const auto depth_comparator = [&](int lhs, int rhs) {
    return ion::math::LengthSquared(point_set.positions[lhs] -
                                    Point3f::Zero()) <
           ion::math::LengthSquared(point_set.positions[rhs] - Point3f::Zero());
  };
  if (partition_count == 1) {
    // In this case, no sorting is actually necessary.
  } else if (partition_count == 2 && point_count >= 3) {
    // This case is very common, and does not require a full sort.  Instead, it
    // is sufficient to ensure that the first and last elements are correct, and
    // that the points are ordered w.r.t. the median.
    //
    // Note that we could easily extend this function to recursively partition
    // point indices via nth_element, but at a certain level, sorting is more
    // efficient.  So, only specialize for the case of two partitions.

    auto median = relevant_point_indices_.begin() + point_count / 2;
    std::nth_element(relevant_point_indices_.begin(), median,
                     relevant_point_indices_.end(), depth_comparator);

    auto minmax = std::minmax_element(relevant_point_indices_.begin(),
                                      relevant_point_indices_.end());
    using std::swap;
    swap(relevant_point_indices_.front(), *minmax.first);
    swap(relevant_point_indices_.back(), *minmax.second);
  } else {
    std::sort(relevant_point_indices_.begin(), relevant_point_indices_.end(),
              depth_comparator);
  }

  for (int i = 0; i < partition_count; ++i) {
    // Reinitialize the partition's GeometryModel based on points which are
    // uniformly-distributed in depth.
    int representative_point = i * point_count / partition_count;
    representative_point = std::min(representative_point, point_count - 1);
    GeometryModel* model = build_partitions[i].GetMutableModel();
    geometry_solver_->InitializeModel(
        relevant_point_indices_[representative_point], model);

    // Assign all points in this depth-range to the partition.

    int low = i * point_count / partition_count;
    int high = (i + 1) * point_count / partition_count;
    high = std::min(high, point_count);
    for (int point_index_index = low; point_index_index < high;
         ++point_index_index) {
      int point_index = relevant_point_indices_[point_index_index];
      float cost = geometry_solver_->ComputeError(point_index, *model);
      build_partitions[i].AddPoint(point_index, cost);
    }
  }

  const int kThreadCount = 1;
  CanonicalizePartitions(kThreadCount, build_partitions);
}

void RobustReinitializingPartitioner::Run(
    const PointSet& point_set, absl::Span<BuildPartition> build_partitions) {
  bool must_reinitialize = false;
  for (const auto& bp : build_partitions) {
    if (!std::isfinite(bp.GetTotalError())) {
      must_reinitialize = true;
      // Ignore first-time "reinitialization".
      if (build_partitions.size() > 1) {
        base::SeuratWarning(
            "Restarting clustering after failure with count = " +
            std::to_string(build_partitions.size()));
      }
      break;
    }
  }
  if (must_reinitialize) {
    reinitializing_stage_->Run(point_set, build_partitions);
  } else {
    regular_stage_->Run(point_set, build_partitions);
  }
}

void HierarchicalPartitioner::Init(const PointSet& point_set) {
  initial_stage_->Init(point_set);
  iterative_stage_->Init(point_set);
}

void HierarchicalPartitioner::Run(const PointSet& point_set,
                                  absl::Span<BuildPartition> build_partitions) {
  const int final_partition_count = build_partitions.size();

  std::vector<int> all_partition_points;
  for (auto& p : build_partitions) {
    std::copy(p.GetPointIndices().begin(), p.GetPointIndices().end(),
              std::back_inserter(all_partition_points));
    p.Clear();
  }

  const int initial_partition_count =
      std::min(initial_partition_count_, final_partition_count);

  // The first N partitions, with N increasing on subsequent iterations.
  absl::Span<BuildPartition> current_partition_slice =
      build_partitions.subspan(0, initial_partition_count);

  // Initialize such that the first partition has all relevant points.
  for (int point_index : all_partition_points) {
    build_partitions.front().AddPoint(point_index, 0);
  }

  initial_stage_->Run(point_set, current_partition_slice);

  while (current_partition_slice.size() < final_partition_count) {
    int partition_count = current_partition_slice.size();
    // Double the number of partitions, but clamp at final_partition_count.
    int num_partitions_to_add =
        std::min(partition_count, final_partition_count - partition_count);
    partition_count += num_partitions_to_add;
    CHECK_GT(num_partitions_to_add, 0);

    current_partition_slice = build_partitions.subspan(0, partition_count);
    iterative_stage_->Run(point_set, current_partition_slice);
  }
  DCHECK_EQ(build_partitions.size(), current_partition_slice.size());
}

}  // namespace tiler
}  // namespace seurat
