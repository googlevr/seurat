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
#include <numeric>
#include <vector>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "absl/random/random.h"
#include "seurat/base/util.h"
#include "seurat/tiler/build_partition.h"
#include "seurat/tiler/geometry_model.h"
#include "seurat/tiler/geometry_solver.h"
#include "seurat/tiler/tiler_test_util.h"

namespace seurat {
namespace tiler {
namespace {

using ion::math::Point3f;
using ion::math::Vector3f;

std::vector<Point3f> BuildTestingPoints(int point_count) {
  std::vector<Point3f> points(point_count);
  for (int i = 0; i < point_count; ++i) {
    points[i] = {static_cast<float>(i), 0.0f, static_cast<float>(i)};
  }
  return points;
}

// Generates a random initial partitioning.
std::vector<BuildPartition> RandomPartitioning(
    const GeometrySolver& geometry_solver, int partition_count,
    int point_count) {
  CHECK_GE(partition_count, 0);
  if (partition_count == 0) {
    return {};
  }

  std::vector<BuildPartition> partitioning(partition_count);

  // Initialize the partitioning randomly by popping point_indices off of a
  // randomly-permuted vector of all indices.
  std::vector<int> all_point_indices(point_count);
  std::iota(all_point_indices.begin(), all_point_indices.end(), 0);
  std::shuffle(all_point_indices.begin(), all_point_indices.end(),
               absl::SharedBitGen());

  std::mt19937 random;
  std::uniform_int_distribution<int> partition_to_assign(
      0, partitioning.size() - 1);
  for (int index = 0; index < point_count; ++index) {
    BuildPartition& p = partitioning[partition_to_assign(random)];
    if (p.Empty()) {
      geometry_solver.InitializeModel(index, p.GetMutableModel());
    }
    p.AddPoint(index, geometry_solver.ComputeError(index, p.GetModel()));
  }
  return partitioning;
}

TEST(PartitionerStagesTest, TestRandomizedPartitionInitializationStage) {
  const int kThreadCount = 1;
  auto geometry_solver = std::make_shared<FakeGeometrySolver>();

  RandomizedPartitionInitializationStage stage(kThreadCount, geometry_solver);

  std::vector<Point3f> points = BuildTestingPoints(1000);
  PointSet point_set = {0, points, {}, {}};

  stage.Init(point_set);

  std::vector<BuildPartition> partitioning(100);
  stage.Run(point_set, absl::MakeSpan(partitioning));

  std::vector<int> initial_seed_surfels;
  for (const BuildPartition& bp : partitioning) {
    EXPECT_TRUE(bp.Empty());
    initial_seed_surfels.push_back(static_cast<int>(bp.GetModel().center[0]));
    EXPECT_EQ(0.0f, bp.GetModel().center[1]);
    EXPECT_EQ(0.0f, bp.GetModel().center[2]);
  }

  std::sort(initial_seed_surfels.begin(), initial_seed_surfels.end());
  auto unique_end =
      std::unique(initial_seed_surfels.begin(), initial_seed_surfels.end());
  int num_unique_seeds = unique_end - initial_seed_surfels.begin();
  EXPECT_EQ(partitioning.size(), num_unique_seeds);
  ExpectNoDuplicatePoints(partitioning);
}

TEST(PartitionerStagesTest, TestSurfaceProxyRefinementStage) {
  const int kThreadCount = 1;
  auto geometry_solver = std::make_shared<FakeGeometrySolver>();

  GeometryModelRefinementStage stage(kThreadCount, geometry_solver);

  std::vector<Point3f> points = BuildTestingPoints(10000);
  PointSet point_set = {0, points, {}, {}};

  stage.Init(point_set);

  std::vector<BuildPartition> partitioning =
      RandomPartitioning(*geometry_solver, 100, points.size());

  float random_partitioning_error =
      ComputeTotalError(*geometry_solver, partitioning);

  stage.Run(point_set, absl::MakeSpan(partitioning));

  float post_refinement_error =
      ComputeTotalError(*geometry_solver, partitioning);
  EXPECT_GE(random_partitioning_error, post_refinement_error);

  ExpectAllPointsPresent(partitioning, points.size());
  ExpectNoDuplicatePoints(partitioning);
}

TEST(PartitionerStagesTest, TestGreedyPointAssignmentStage) {
  const int kThreadCount = 1;
  auto geometry_solver = std::make_shared<FakeGeometrySolver>();

  GreedyPointAssignmentStage stage(kThreadCount, 16, geometry_solver);

  std::vector<Point3f> points = BuildTestingPoints(10000);
  PointSet point_set = {0, points, {}, {}};

  stage.Init(point_set);

  std::vector<BuildPartition> partitioning =
      RandomPartitioning(*geometry_solver, 100, points.size());

  float random_partitioning_error =
      ComputeTotalError(*geometry_solver, partitioning);

  stage.Run(point_set, absl::MakeSpan(partitioning));

  float post_refinement_error =
      ComputeTotalError(*geometry_solver, partitioning);
  EXPECT_GE(random_partitioning_error, post_refinement_error);

  ExpectAllPointsPresent(partitioning, points.size());
  ExpectNoDuplicatePoints(partitioning);
}

TEST(PartitionerStagesTest,
     TestGreedyPointAssignmentStage_InvalidProxyFitting) {
  // Test point-partition assignment when the GeometrySolver cost function
  // returns infinite error for all pairings.
  //
  // The precise behavior of how to handle this edge-case is subject to change
  // (this test currently asserts that such points are omitted from the
  // partitioning).
  //
  // However, at the very least, it should not crash.
  const int kThreadCount = 1;
  auto geometry_solver = std::make_shared<InvalidGeometrySolver>();

  std::vector<Point3f> points = BuildTestingPoints(10000);
  PointSet point_set = {0, points, {}, {}};

  FakeGeometrySolver valid_cost_function;
  valid_cost_function.Init(point_set);
  std::vector<BuildPartition> partitioning =
      RandomPartitioning(valid_cost_function, 100, points.size());

  GreedyPointAssignmentStage stage(kThreadCount, 16, geometry_solver);
  stage.Init(point_set);
  stage.Run(point_set, absl::MakeSpan(partitioning));

  for (const auto& bp : partitioning) {
    EXPECT_TRUE(bp.Empty());
  }
}

TEST(PartitionerStagesTest, TestPartitionSplittingStage) {
  auto geometry_solver = std::make_shared<FakeGeometrySolver>();

  std::vector<Point3f> points = BuildTestingPoints(10000);
  PointSet point_set = {0, points, {}, {}};

  PartitionSplittingStage stage(geometry_solver);
  stage.Init(point_set);

  std::vector<BuildPartition> original_partitioning =
      RandomPartitioning(*geometry_solver, 100, points.size());
  std::vector<BuildPartition> partitioning = original_partitioning;

  // Run the stage and verify that the resulting partitioning is as expected.
  stage.Run(point_set, absl::MakeSpan(partitioning));
  ExpectAllPointsPresent(partitioning, points.size());
  EXPECT_EQ(100, partitioning.size());
  EXPECT_GE(ComputeTotalError(*geometry_solver, original_partitioning),
            ComputeTotalError(*geometry_solver, partitioning));
  for (const auto& bp : partitioning) {
    EXPECT_FALSE(bp.Empty());
  }

  // Running a second time should do nothing.
  std::vector<BuildPartition> partitioning2 = partitioning;
  stage.Run(point_set, absl::MakeSpan(partitioning));
  EXPECT_EQ(partitioning2, partitioning);

  // Add more (empty) partitions and run again.
  for (int i = 0; i < 20; ++i) {
    partitioning.push_back(BuildPartition{});
  }
  stage.Run(point_set, absl::MakeSpan(partitioning));
  ExpectAllPointsPresent(partitioning, points.size());
  ExpectNoDuplicatePoints(partitioning);
  EXPECT_EQ(100 + 20, partitioning.size());
  EXPECT_GE(ComputeTotalError(*geometry_solver, original_partitioning),
            ComputeTotalError(*geometry_solver, partitioning));
  for (const auto& bp : partitioning) {
    EXPECT_FALSE(bp.Empty());
  }
}

TEST(PartitionerStagesTest, PointExchangeStage) {
  const int kThreadCount = 1;
  auto geometry_solver = std::make_shared<FakeGeometrySolver>();

  std::vector<Point3f> points = BuildTestingPoints(10000);
  PointSet point_set = {0, points, {}, {}};

  PointExchangeStage stage(kThreadCount, geometry_solver);
  stage.Init(point_set);

  const int kNumPointsOfInterest = 500;
  const int kPartitionCount = 15;

  // Construct a random partitioning which only includes the first
  // kNumPointsOfInterest points.
  std::vector<BuildPartition> original_partitioning = RandomPartitioning(
      *geometry_solver, kPartitionCount, kNumPointsOfInterest);
  std::vector<BuildPartition> partitioning = original_partitioning;

  // Run the stage and verify that the resulting partitioning is as expected.
  stage.Run(point_set, absl::MakeSpan(partitioning));
  ExpectAllPointsPresent(partitioning, kNumPointsOfInterest);
  EXPECT_EQ(kPartitionCount, partitioning.size());
  for (const auto& bp : partitioning) {
    EXPECT_FALSE(bp.Empty());
  }
  EXPECT_GE(ComputeTotalError(*geometry_solver, original_partitioning),
            ComputeTotalError(*geometry_solver, partitioning));

  // Running a second time should do nothing.
  std::vector<BuildPartition> partitioning2 = partitioning;
  stage.Run(point_set, absl::MakeSpan(partitioning));
  EXPECT_EQ(partitioning2, partitioning);
}

void TestDepthBasedRedistributionStage(int initial_partition_count) {
  auto geometry_solver = std::make_shared<FakeGeometrySolver>();

  std::vector<Point3f> points = BuildTestingPoints(10000);
  PointSet point_set = {0, points, {}, {}};

  DepthBasedRedistributionStage stage(geometry_solver);
  stage.Init(point_set);

  const int kNumPointsOfInterest = 500;

  // Construct a random partitioning which only includes the first
  // kNumPointsOfInterest points.
  std::vector<BuildPartition> original_partitioning = RandomPartitioning(
      *geometry_solver, initial_partition_count, kNumPointsOfInterest);
  std::vector<BuildPartition> partitioning = original_partitioning;

  // Run the stage and verify that the resulting partitioning is as expected.
  stage.Run(point_set, absl::MakeSpan(partitioning));
  ExpectAllPointsPresent(partitioning, kNumPointsOfInterest);
  EXPECT_EQ(initial_partition_count, partitioning.size());
  for (const auto& bp : partitioning) {
    EXPECT_FALSE(bp.Empty());
  }

  // Running a second time should do nothing.
  std::vector<BuildPartition> partitioning2 = partitioning;
  stage.Run(point_set, absl::MakeSpan(partitioning));
  EXPECT_EQ(partitioning2, partitioning);

  // Add more (empty) partitions and run again.
  for (int i = 0; i < 20; ++i) {
    partitioning.push_back(BuildPartition{});
  }
  stage.Run(point_set, absl::MakeSpan(partitioning));
  ExpectAllPointsPresent(partitioning, kNumPointsOfInterest);
  ExpectNoDuplicatePoints(partitioning);
  EXPECT_EQ(initial_partition_count + 20, partitioning.size());
  for (const auto& bp : partitioning) {
    EXPECT_FALSE(bp.Empty());
  }
}

TEST(PartitionerStagesTest, TestDepthBasedRedistributionStage_1_Partition) {
  TestDepthBasedRedistributionStage(1);
}

TEST(PartitionerStagesTest, TestDepthBasedRedistributionStage_2_Partitions) {
  TestDepthBasedRedistributionStage(2);
}

TEST(PartitionerStagesTest, TestDepthBasedRedistributionStage_15_Partitions) {
  TestDepthBasedRedistributionStage(15);
}

TEST(PartitionerStagesTest, TestRobustReinitializingPartitioner) {
  auto fake_reinitializer_stage = std::make_shared<TestingPartitionerStage>();
  auto fake_regular_stage = std::make_shared<TestingPartitionerStage>();
  RobustReinitializingPartitioner robust_partitioner(fake_reinitializer_stage,
                                                     fake_regular_stage);

  std::vector<Point3f> points = BuildTestingPoints(10000);
  PointSet point_set = {0, points, {}, {}};
  FakeGeometrySolver geometry_solver;
  geometry_solver.Init(point_set);
  std::vector<BuildPartition> partitioning =
      RandomPartitioning(geometry_solver, 100, points.size());

  robust_partitioner.Init(point_set);
  EXPECT_EQ(1, fake_regular_stage->GetInitCount());
  robust_partitioner.Run(point_set, absl::MakeSpan(partitioning));
  EXPECT_EQ(1, fake_regular_stage->GetRunCount());
  EXPECT_EQ(0, fake_reinitializer_stage->GetRunCount());

  partitioning[0].AddPoint(0, std::numeric_limits<float>::infinity());
  robust_partitioner.Run(point_set, absl::MakeSpan(partitioning));
  EXPECT_EQ(1, fake_regular_stage->GetRunCount());
  EXPECT_EQ(1, fake_reinitializer_stage->GetRunCount());
}

TEST(PartitionerStagesTest, TestCompositeStages) {
  std::vector<Point3f> points = BuildTestingPoints(10000);
  PointSet point_set = {0, points, {}, {}};

  std::vector<BuildPartition> partitioning(100);

  // Build up a composition of multiple stages as follows:
  // SequentialPartitioner:
  //  1. TestingPartitionStage #1
  //  2. IterativePartitioner:
  //       * TestingPartitionStage #1
  //       * TestingPartitionStage #1
  //       * TestingPartitionStage #1
  //  3. TestingPartitionStage #2
  auto testing_stage1 = std::make_shared<TestingPartitionerStage>();
  auto testing_stage2 = std::make_shared<TestingPartitionerStage>();

  auto repeated1 = std::make_shared<IterativePartitioner>(3, testing_stage1);
  auto sequence = std::make_shared<SequentialPartitioner>(
      std::vector<std::shared_ptr<PartitionerStage>>{testing_stage1, repeated1,
                                                     testing_stage2});

  sequence->Init(point_set);

  EXPECT_LE(1, testing_stage1->GetInitCount());
  EXPECT_GE(2, testing_stage1->GetInitCount());
  EXPECT_EQ(1, testing_stage2->GetInitCount());

  sequence->Run(point_set, absl::MakeSpan(partitioning));

  EXPECT_EQ(4, testing_stage1->GetRunCount());
  EXPECT_EQ(1, testing_stage2->GetRunCount());
}

TEST(PartitionerStagesTest, TestHierarchicalPartitioner) {
  auto geometry_solver = std::make_shared<FakeGeometrySolver>();
  const int kThreadCount = 1;

  auto random_initialization =
      std::make_shared<RandomizedPartitionInitializationStage>(kThreadCount,
                                                               geometry_solver);
  auto splitting = std::make_shared<PartitionSplittingStage>(geometry_solver);
  std::shared_ptr<PartitionerStage> fast_assignment(
      new GreedyPointAssignmentStage(kThreadCount, 16, geometry_solver));
  auto refitting = std::make_shared<GeometryModelRefinementStage>(
      kThreadCount, geometry_solver);

  std::shared_ptr<PartitionerStage> initial_stage(
      new SequentialPartitioner({random_initialization,  //
                                 fast_assignment,        //
                                 refitting}));
  std::shared_ptr<PartitionerStage> iterated_stage(
      new SequentialPartitioner({splitting,        //
                                 fast_assignment,  //
                                 refitting}));

  const int kInitialPartitionCount = 100;
  const int kFinalPartitionCount = 100;
  auto partitioner = std::make_shared<HierarchicalPartitioner>(
      kInitialPartitionCount, initial_stage, iterated_stage);

  std::vector<Point3f> points = BuildTestingPoints(1000);
  PointSet point_set = {0, points, {}, {}};

  std::vector<BuildPartition> partitioning(kFinalPartitionCount,
                                           BuildPartition{});
  partitioner->Init(point_set);
  partitioner->Run(point_set, absl::MakeSpan(partitioning));

  ExpectAllPointsPresent(partitioning, points.size());
  ExpectNoDuplicatePoints(partitioning);
  EXPECT_EQ(kFinalPartitionCount, partitioning.size());
}

}  // namespace
}  // namespace tiler
}  // namespace seurat
