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

#include "seurat/tiler/tiler.h"

#include "ion/math/range.h"
#include "seurat/tiler/build_partition.h"
#include "seurat/tiler/candidate_tile_generator.h"
#include "seurat/tiler/partitioner_stage.h"
#include "seurat/tiler/partitioner_stages.h"
#include "seurat/tiler/rail_disk_solver.h"
#include "seurat/tiler/selection/solver_factory.h"
#include "seurat/tiler/selection_tiler.h"
#include "seurat/tiler/subdivision.h"
#include "seurat/tiler/subset_geometry_solver.h"

namespace seurat {
namespace tiler {

using ion::math::Range1f;

namespace {

std::unique_ptr<CandidateTileGenerator> MakeCandidateTileGenerator(
    const TilerFactory::Parameters& parameters,
    const std::shared_ptr<Subdivision>& subdivision) {
  // The number of points to select within each bin when fitting planes.
  //
  // If INT_MAX, this uses all points.
  //
  // Decreasing this increases speed, at a loss of quality.
  const int kMaxPointsPerDiskFittingOptimization = 250;

  // The maximum number of partitions per bin.
  //
  // Higher allows for more flexibility in the solver, but beyond 8 is likely
  // past the point of diminishing returns.
  //
  // Running-time is roughly-linear in this parameter.
  const int kMaxPartitionsPerBin = 8;

  // The range of tile corner depths (distances from the origin) to permit
  // without penalty.
  //
  // Note that this is *not* a hard limit on the depth-range of
  // resulting tiles.
  //
  // Scale the near limit (by 0.01) to support grazing-angle geometry which may
  // enter the headbox.
  //
  // Scale the far limit by sqrt(3), since this depth_range is applied uniformly
  // to ray-distances, as opposed to using an explicit origin-centered AABB
  // model.
  const float skybox_bounding_sphere_radius =
      std::sqrt(3.0f) * parameters.skybox_radius;
  const Range1f depth_range(0.01f * parameters.headbox_radius,
                            skybox_bounding_sphere_radius);

  // Early iterations of the partitioner will include the tangential-term,
  // resulting in a more-convex problem.
  //
  // The final iteration will ignore the tangential-term, since it better
  // reflects the visual distortion of the final geometry.
  //
  // This is analogous to the idea of graduated-nonconvexity.
  auto disk_solver = [=]() -> std::shared_ptr<GeometrySolver> {
    const float kTangentialFactor = 0.01f;
    return std::make_shared<RailDiskSolver>(kTangentialFactor, subdivision,
                                            depth_range);
  };
  auto plane_solver = [=]() -> std::shared_ptr<GeometrySolver> {
    const float kTangentialFactor = 0.0f;
    return std::make_shared<RailDiskSolver>(kTangentialFactor, subdivision,
                                            depth_range);
  };

  auto approximate_disk_solver = [=]() {
    return std::make_shared<SubsetGeometrySolver>(
        kMaxPointsPerDiskFittingOptimization, disk_solver());
  };

  auto tile_resolver = [=]() -> std::shared_ptr<TileResolver> {
    return std::make_shared<RailTileResolver>(subdivision);
  };

  auto candidate_generator = [=]() -> std::shared_ptr<CandidateTileGenerator> {
    auto initialization = std::make_shared<DepthBasedRedistributionStage>(
        approximate_disk_solver());
    auto splitting =
        std::make_shared<PartitionSplittingStage>(approximate_disk_solver());
    auto assignment =
        std::make_shared<PointExchangeStage>(1, approximate_disk_solver());
    auto refitting = std::make_shared<GeometryModelRefinementStage>(
        1, approximate_disk_solver());
    auto final_assignment =
        std::make_shared<PointExchangeStage>(1, plane_solver());

    // The regular stage performs a sequence of assignment, splitting, refitting
    // steps.
    //
    // Note that all steps use the disk_fitting objective-function, except the
    // final_assignment stage which uses the disk_solver objective.
    //
    // This ensures that the candidate partitions do not include the "tangential
    // term" when sent to the PyramidPartitioner (see below) for selection.  So,
    // the overall Partitioner optimizes for the disk_solver model, despite
    // inner loop using the disk_fitting metric.
    std::shared_ptr<PartitionerStage> regular_stage(
        new SequentialPartitioner({splitting,   //
                                   refitting,   //
                                   assignment,  //
                                   refitting,   //
                                   final_assignment}));

    std::shared_ptr<PartitionerStage> initial_stage(
        new SequentialPartitioner({initialization,  //
                                   refitting}));
    std::shared_ptr<PartitionerStage> reinitializing_stage =
        std::make_shared<HierarchicalPartitioner>(
            2, initial_stage,
            std::shared_ptr<PartitionerStage>(
                new SequentialPartitioner({splitting,   //
                                           refitting,   //
                                           assignment,  //
                                           refitting,   //
                                           final_assignment})));

    // To ensure robustness when generating candidate partitionings, use the
    // regular_stage by default, but switch to reinitializing_stage if we run
    // into numerical stability problems.
    auto robustifier = std::make_shared<RobustReinitializingPartitioner>(
        reinitializing_stage, regular_stage);

    return std::make_shared<ExhaustiveCandidateTileGenerator>(
        kMaxPartitionsPerBin, robustifier, plane_solver(), tile_resolver());
  };

  return std::unique_ptr<CandidateTileGenerator>(
      new ParallelCandidateTileGenerator(parameters.thread_count,
                                         candidate_generator));
}

}  // namespace

// static
std::unique_ptr<Tiler> TilerFactory::CreateSelectionTiler(
    const Parameters& parameters) {
  // High level flow:
  //
  // SelectionTiler:
  //  * Generate candidates using the CandidateTileGenerator:
  //      * Invoke a (complex) sequence of PartitionerStages to refine
  //        BuildPartitions.
  //      * Resolve the BuildPartitions into Tiles via the TileResolver.
  //  * Generate a SelectionProblem encoding the selection of these tiles, using
  //    weights computed from the TileWeightModel to constrain the selection.
  //  * Invoke the SelectionSolver.
  //  * Return the selected tiles.

  // TODO(puneetl):  Don't hardcode dilation parameters.
  const int kInputPixelsPerDegree = 11;
  const float kDilationFactor = 1.5f;
  const float kDilationRadians =
      kDilationFactor * 2.0f * M_PI / (kInputPixelsPerDegree * 360.0f);

  auto subdivision = std::make_shared<BoundsDilatingSubdivision>(
      kDilationRadians,
      std::unique_ptr<Subdivision>(
          new CubemapQuadtreeSubdivision(parameters.max_subdivision_level)));

  selection::SelectionSolverParameters selection_params;
  selection_params.thread_count = parameters.thread_count;
  std::shared_ptr<selection::SelectionSolver> selection_solver(
      selection::CreateSelectionSolver(selection_params));

  std::vector<std::unique_ptr<TileWeightModel>> weight_models;
  weight_models.emplace_back(new TriangleCountTileWeightModel);
  weight_models.emplace_back(new ProjectedAreaTileWeightModel);
  weight_models.emplace_back(DirectionalOverdrawTileWeightModel::Build(
      parameters.peak_overdraw_samples,
      parameters.peak_overdraw_field_of_view_degrees * (2.0f * M_PI) / 360.0f,
      parameters.headbox_radius));
  auto tile_weight_model =
      std::make_shared<CombinedTileWeightModel>(std::move(weight_models));

  // This vector must match weight_models.
  std::vector<float> max_weight;
  max_weight.push_back(parameters.tile_count * 2.0f);
  max_weight.push_back(parameters.overdraw_factor);
  for (int i = 0; i < parameters.peak_overdraw_samples; ++i) {
    max_weight.push_back(parameters.peak_overdraw_factor);
  }

  return std::unique_ptr<Tiler>(new SelectionTiler(
      parameters.min_subdivision_level, subdivision,
      std::shared_ptr<CandidateTileGenerator>(
          MakeCandidateTileGenerator(parameters, subdivision).release()),
      selection_solver, tile_weight_model, max_weight,
      parameters.thread_count));
}

}  // namespace tiler
}  // namespace seurat
