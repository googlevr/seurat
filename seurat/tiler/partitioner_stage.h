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

#ifndef VR_SEURAT_TILER_PARTITIONER_STAGE_H_
#define VR_SEURAT_TILER_PARTITIONER_STAGE_H_

#include <functional>
#include <vector>

#include "absl/types/span.h"
#include "seurat/tiler/build_partition.h"
#include "seurat/tiler/point_set.h"

namespace seurat {
namespace tiler {

// A partitioning problem is the task of partitioning a set of points into
// disjoint subsets (partitions), each with parameters for a GeometryModel used
// to later render those points.
//
// A PartitionerStage may be composed with other stages to form a solver for a
// partitioning problem which minimizes a particular error metric, e.g. by
// alternating between assigning points to partitions and optimizing each
// partition's GeometryModel.
//
// In general, implementations are not thread-safe, and may internally cache
// intermediate representations, acceleration structures, and scratch space.
class PartitionerStage {
 public:
  // Constructs a new instance of a PartitionerStage.
  using Factory = std::function<std::shared_ptr<PartitionerStage>(void)>;

  virtual ~PartitionerStage() = default;

  // This must always be called before Run().  Partitioner stages may use this
  // to precompute acceleration-structures based on the points to be
  // partitioned.
  //
  // For the sake of efficiency, a single PartitionerStage may be reused to
  // solve multiple different partitioning problems, e.g. to cache allocations
  // for internal data-structures.
  virtual void Init(const PointSet& point_set) {
    // Do nothing by default.
  }

  // Runs the partitioning stage, modifying the |build_partitions|.  The
  // |point_set| is provided again, and must be the same as that previously
  // passed to Init().
  virtual void Run(const PointSet& point_set,
                   absl::Span<BuildPartition> build_partitions) = 0;
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_PARTITIONER_STAGE_H_
