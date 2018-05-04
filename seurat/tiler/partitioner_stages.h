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

#ifndef VR_SEURAT_TILER_PARTITIONER_STAGES_H_
#define VR_SEURAT_TILER_PARTITIONER_STAGES_H_

#include <memory>
#include <random>
#include <vector>

#include "ion/math/vector.h"
#include "seurat/tiler/build_partition.h"
#include "seurat/tiler/geometry_solver.h"
#include "seurat/tiler/partitioner_stage.h"
#include "seurat/tiler/point_set.h"

namespace seurat {
namespace tiler {

// Clears and initializes all partitions with parameters derived from
// randomly-selected points.
//
// After this stage, all partitions will be empty and have new GeometryModel
// parameters.
//
// This stage considers all points in the PointSet.
class RandomizedPartitionInitializationStage : public PartitionerStage {
 public:
  RandomizedPartitionInitializationStage(
      int thread_count, std::shared_ptr<GeometrySolver> geometry_solver)
      : thread_count_(thread_count),
        geometry_solver_(std::move(geometry_solver)) {}
  ~RandomizedPartitionInitializationStage() override = default;

  void Init(const PointSet& point_set) override {
    geometry_solver_->Init(point_set);
  }
  void Run(const PointSet& points,
           absl::Span<BuildPartition> build_partitions) override;

 private:
  const int thread_count_;

  std::mt19937 random_;

  std::shared_ptr<GeometrySolver> geometry_solver_;
};

// Refits the GeometryModel of each partition based on the points currently
// assigned to it.
//
// This stage only considers points in the existing partitioning.
class GeometryModelRefinementStage : public PartitionerStage {
 public:
  GeometryModelRefinementStage(int thread_count,
                               std::shared_ptr<GeometrySolver> geometry_solver)
      : thread_count_(thread_count),
        geometry_solver_(std::move(geometry_solver)) {}
  ~GeometryModelRefinementStage() override = default;

  void Init(const PointSet& point_set) override {
    geometry_solver_->Init(point_set);
  }
  void Run(const PointSet& points,
           absl::Span<BuildPartition> build_partitions) override;

 private:
  const int thread_count_;

  std::shared_ptr<GeometrySolver> geometry_solver_;
};

// For each point, finds neighboring partitions and assigns the point to the
// partition with the lowest cost.
//
// Neighboring partitions are those partitions which are closest to the point
// w.r.t. their direction vectors from the origin.
//
// This stage considers all points in the PointSet.
class GreedyPointAssignmentStage : public PartitionerStage {
 public:
  GreedyPointAssignmentStage(int thread_count, int neighboring_partition_count,
                             std::shared_ptr<GeometrySolver> geometry_solver)
      : thread_count_(thread_count),
        neighboring_partition_count_(neighboring_partition_count),
        geometry_solver_(std::move(geometry_solver)) {}
  ~GreedyPointAssignmentStage() override = default;

  void Init(const PointSet& point_set) override {
    geometry_solver_->Init(point_set);
  }
  void Run(const PointSet& point_set,
           absl::Span<BuildPartition> build_partitions) override;

 private:
  const int thread_count_;
  const int neighboring_partition_count_;
  std::shared_ptr<GeometrySolver> geometry_solver_;
};

// Finds all empty partitions and uses them to subdivide non-empty partitions
// with highest total point-partition cost.
//
// If there are no empty partitions, this stage does nothing.
//
// This stage only considers points in the existing partitioning.
class PartitionSplittingStage : public PartitionerStage {
 public:
  explicit PartitionSplittingStage(
      std::shared_ptr<GeometrySolver> geometry_solver)
      : minimum_points_per_partition_(2),
        geometry_solver_(std::move(geometry_solver)) {}
  ~PartitionSplittingStage() override = default;

  void Init(const PointSet& point_set) override {
    geometry_solver_->Init(point_set);
  }
  void Run(const PointSet& point_set,
           absl::Span<BuildPartition> build_partitions) override;

 private:
  // The minimum size of a partition to be split.  Partitions with fewer points
  // than this are never subdivided.
  const int minimum_points_per_partition_;

  std::shared_ptr<GeometrySolver> geometry_solver_;
};

// Considers only the points in the current partitioning, re-assigning them to
// the partitions with lowest cost.
//
// This stage only considers points in the existing partitioning.
class PointExchangeStage : public PartitionerStage {
 public:
  PointExchangeStage(int thread_count,
                     std::shared_ptr<GeometrySolver> geometry_solver)
      : thread_count_(thread_count),
        geometry_solver_(std::move(geometry_solver)) {}
  ~PointExchangeStage() override = default;

  void Init(const PointSet& point_set) override {
    geometry_solver_->Init(point_set);
  }
  void Run(const PointSet& point_set,
           absl::Span<BuildPartition> build_partitions) override;

 private:
  // Find the partition, from all partitions, which has the least error,
  // returning its index.
  //
  // Use the first partition by default, if all partitions have infinite
  // error.  In practice, this rarely happens, and other stages (e.g.
  // RobustReinitializationPartitioner) can compensate for this.
  int FindBestPartition(absl::Span<const BuildPartition> partitions,
                        int point_index, float* cost) const;

  const int thread_count_;
  const std::shared_ptr<GeometrySolver> geometry_solver_;

  // A temporary vector storing the indices of all points to redistribute.
  //
  // Cached here to reuse temporary allocations.
  std::vector<int> relevant_point_indices_;
};

// Like PointExchangeStage, but (re)initializes all partitions to have geometry
// models generated from equally-spaced points (based on depth ordering).
//
// In other words, this considers all points in the current partitioning, sorts
// them based on depth, and splits those points into equal-sized partitions.
//
// The end result is that all points in all partitions are redistributed based
// on their depth, and new geometry models are initialized.
//
// This stage only considers points in the existing partitioning.
class DepthBasedRedistributionStage : public PartitionerStage {
 public:
  explicit DepthBasedRedistributionStage(
      std::shared_ptr<GeometrySolver> geometry_solver)
      : geometry_solver_(std::move(geometry_solver)) {}
  ~DepthBasedRedistributionStage() override = default;

  void Init(const PointSet& point_set) override {
    geometry_solver_->Init(point_set);
  }
  void Run(const PointSet& point_set,
           absl::Span<BuildPartition> build_partitions) override;

 private:
  const std::shared_ptr<GeometrySolver> geometry_solver_;

  // A temporary vector storing the indices of all points to redistribute.
  //
  // Cached here to reuse temporary allocations.
  std::vector<int> relevant_point_indices_;
};

// Processes the existing partitioning as follows:
//
// If the partitioning is ill-conditioned (any of the partitions have infinite
// cost), then the entire partitioning is reinitialized using the
// reinitialization_stage.
//
// Otherwise, the regular_stage is run.
class RobustReinitializingPartitioner : public PartitionerStage {
 public:
  RobustReinitializingPartitioner(
      std::shared_ptr<PartitionerStage> reinitializing_stage,
      std::shared_ptr<PartitionerStage> regular_stage)
      : reinitializing_stage_(std::move(reinitializing_stage)),
        regular_stage_(std::move(regular_stage)) {}

  void Init(const PointSet& point_set) override {
    reinitializing_stage_->Init(point_set);
    regular_stage_->Init(point_set);
  }
  void Run(const PointSet& point_set,
           absl::Span<BuildPartition> build_partitions) override;

 private:
  // The stage to run if the current partitioning has numerical stability
  // problems (has a partition with non-finite cost).
  std::shared_ptr<PartitionerStage> reinitializing_stage_;

  // The stage to run if the current partitioning is okay.
  std::shared_ptr<PartitionerStage> regular_stage_;
};

// Executes a sequence of PartitionerStages in order.
class SequentialPartitioner : public PartitionerStage {
 public:
  explicit SequentialPartitioner(
      std::vector<std::shared_ptr<PartitionerStage>> child_stages)
      : child_stages_(std::move(child_stages)) {}

  ~SequentialPartitioner() override = default;

  void Init(const PointSet& point_set) override {
    for (auto& child : child_stages_) child->Init(point_set);
  }
  void Run(const PointSet& point_set,
           absl::Span<BuildPartition> build_partitions) override {
    for (auto& child : child_stages_) child->Run(point_set, build_partitions);
  }

 private:
  std::vector<std::shared_ptr<PartitionerStage>> child_stages_;
};

// Repeatedly executes a PartitionerStage for a predetermined number of
// iterations.
class IterativePartitioner : public PartitionerStage {
 public:
  IterativePartitioner(int iteration_count,
                       std::shared_ptr<PartitionerStage> child)
      : iteration_count_(iteration_count), child_(std::move(child)) {}

  ~IterativePartitioner() override = default;

  void Init(const PointSet& point_set) override { child_->Init(point_set); }
  void Run(const PointSet& point_set,
           absl::Span<BuildPartition> build_partitions) override {
    for (int i = 0; i < iteration_count_; ++i) {
      child_->Run(point_set, build_partitions);
    }
  }

 private:
  const int iteration_count_;

  const std::shared_ptr<PartitionerStage> child_;
};

// A partitioner which executes the specified other PartitionerStages,
// doubling the number of partitions each time until all partitions are used.
//
// |initial_stage| is run the first time, and |iterative_stage| is executed
// for subsequent iterations, until all partitions are used.
class HierarchicalPartitioner : public PartitionerStage {
 public:
  HierarchicalPartitioner(int initial_partition_count,
                          std::shared_ptr<PartitionerStage> initial_stage,
                          std::shared_ptr<PartitionerStage> iterative_stage)
      : initial_partition_count_(initial_partition_count),
        initial_stage_(std::move(initial_stage)),
        iterative_stage_(std::move(iterative_stage)) {
    CHECK_GE(initial_partition_count, 0);
  }
  ~HierarchicalPartitioner() override = default;

  void Init(const PointSet& point_set) override;
  void Run(const PointSet& point_set,
           absl::Span<BuildPartition> build_partitions) override;

 private:
  const int initial_partition_count_;
  std::shared_ptr<PartitionerStage> initial_stage_;
  std::shared_ptr<PartitionerStage> iterative_stage_;
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_PARTITIONER_STAGES_H_
