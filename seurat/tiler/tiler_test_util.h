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

#ifndef VR_SEURAT_TILER_TILER_TEST_UTIL_H_
#define VR_SEURAT_TILER_TILER_TEST_UTIL_H_

#include <vector>

#include "ion/math/vector.h"
#include "seurat/tiler/build_partition.h"
#include "seurat/tiler/geometry_model.h"
#include "seurat/tiler/geometry_solver.h"
#include "seurat/tiler/partitioner_stage.h"
#include "seurat/tiler/point_set.h"
#include "seurat/tiler/subdivision.h"
#include "seurat/tiler/tiler.h"
#include "seurat/tiler/tiler_test_util.h"

namespace seurat {
namespace tiler {

// A very simple GeometrySolver implementation which optimizes GeometryModel
// instances such that the absolute distance of each point, according to it's
// *index* to the x-coordinate of the GeometryModel center is minimized.
class FakeGeometrySolver : public GeometrySolver {
 public:
  ~FakeGeometrySolver() override = default;

  // GeometrySolver implementation.
  void Init(const PointSet& point_set) override;
  void InitializeModel(int point_index, GeometryModel* model) const override;
  bool FitModel(absl::Span<const int> point_indices,
                GeometryModel* model) const override;
  float ComputeError(int point_index,
                     const GeometryModel& model) const override;

 private:
  PointSet point_set_;
};

// A GeometrySolver implementation which returns infinite error.
class InvalidGeometrySolver : public FakeGeometrySolver {
 public:
  float ComputeError(int point_index,
                     const GeometryModel& model) const override;
};

// A PartitionerStage which tracks how many times it has been invoked.
class TestingPartitionerStage : public PartitionerStage {
 public:
  ~TestingPartitionerStage() override = default;

  int GetInitCount() const { return init_count_; }
  int GetRunCount() const { return run_count_; }
  PointSet GetPointSet() const { return point_set_; }

  // PartitionerStage implementation.
  void Init(const PointSet& point_set) override;
  void Run(const PointSet& point_set,
           absl::Span<BuildPartition> build_partitions) override;

 private:
  PointSet point_set_;
  int init_count_;
  int run_count_;
};

// Asserts that all points are assigned to exactly one partition.
void ExpectAllPointsPresent(const std::vector<BuildPartition>& partitioning,
                            int point_count);

// Asserts that there are no duplicate points in the given BuildPartitions.
void ExpectNoDuplicatePoints(const std::vector<BuildPartition>& partitioning);

// Computes the total error of a given |partitioning|.
float ComputeTotalError(const GeometrySolver& solver,
                        const std::vector<BuildPartition>& partitioning);

// Returns |num_points| randomly-selected points on the unit-sphere.
std::vector<ion::math::Point3f> GenerateUnitSpherePoints(int num_points);

// Returns true if rays from the origin to all points intersect a tile.
void ExpectTilesCoverAllPoints(absl::Span<const Tile> tiles,
                               absl::Span<const ion::math::Point3f> points);

// Uses the given |tiler| to partition a set of points on a spherical cap
// consisting of all points, p, on the unit sphere such that
//   (p DOT direction) > min_dot
//
// The resulting tiles are expected to be valid:  Rays from the origin to all
// points must intersect a tile.
void ExpectTilesSphericalCap(Tiler* tiler, int point_set_id, int point_count,
                             const ion::math::Vector3f& direction,
                             float min_dot);

// Uses the given |partitioner| to partition a set of points on the unit sphere.
//
// The resulting tiles are expected to be valid:  Rays from the origin to all
// points must intersect a tile.
void ExpectTilesUnitSphere(Tiler* tiler, int point_set_id, int point_count);

// Validates the given |subdivision|, ensuring it contains all points and that
// each child node contains only a subset of its parent's points.
void ExpectWellFormedSubdivision(const PointSet& points,
                                 const Subdivision& subdivision, int max_depth);

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_TILER_TEST_UTIL_H_
