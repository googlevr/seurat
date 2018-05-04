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

// Framework for evaluating Seurat Artifacts (e.g. diff'ing rendered images,
// evaluating simulated runtime rendering cost).
//
// * ArtifactEvaluator provides an interface for plugging in different
//   evaluators.
// * EvaluationExporter packages these evaluators as an Exporter
//   which "exports" the human-readable results (writes them to a file).
//
// Various evaluators are implemented:
// * OverdrawEvaluator - Wraps a CostEstimator to evaluate the overdraw required
//   to render a mesh.
// * RenderEvaluator - Wraps a RenderSim to evaluate visual quality relative to
//   the original images.
#ifndef VR_SEURAT_ARTIFACT_EVALUATION_ARTIFACT_EVALUATOR_H_
#define VR_SEURAT_ARTIFACT_EVALUATION_ARTIFACT_EVALUATOR_H_

#include <memory>

#include "ion/math/fieldofview.h"
#include "ion/math/vector.h"
#include "seurat/artifact/artifact.h"
#include "seurat/artifact/artifact_processor.h"
#include "seurat/artifact/evaluation/cost_estimator.h"
#include "seurat/artifact/evaluation/render_sim.h"
#include "seurat/base/color.h"
#include "seurat/base/file_system.h"
#include "seurat/image/image.h"
#include "seurat/ingest/view_group_loader.h"

namespace seurat {
namespace artifact {

// An interface for methods of evaluating quality of an Artifact.
//
// For example, an implementation may render out the Artifact and compare the
// results to the original views, resulting in summary statistics of
// image-comparison metrics.
class ArtifactEvaluator {
 public:
  // A human-readable name & value for a particular metric which was evaluated.
  struct Result {
    std::string name;
    std::string value;
  };

  virtual ~ArtifactEvaluator() = default;

  // The human-readable name of the evaluator.
  virtual std::string Name() const = 0;

  // Evaluates the quality of the |artifact|.
  virtual std::vector<Result> Evaluate(
      const artifact::Artifact& artifact) const = 0;
};

// Evaluates an Artifact and exports the results to file.
//
// This is a simple adaptor to make an ArtifactEvaluator act like an Exporter.
class EvaluationExporter : public ArtifactProcessor {
 public:
  EvaluationExporter(std::shared_ptr<base::FileSystem> filesystem,
                     std::string basename,
                     std::vector<std::unique_ptr<ArtifactEvaluator>> evaluators)
      : filesystem_(std::move(filesystem)),
        basename_(std::move(basename)),
        evaluators_(std::move(evaluators)) {}
  ~EvaluationExporter() override = default;

  // Exporter implementation.
  base::Status Process(artifact::Artifact* artifact) const override;

 private:
  // Used to write output files.
  const std::shared_ptr<base::FileSystem> filesystem_;

  // The basename for the .eval text file to write.
  const std::string basename_;

  // The evaluators to run.
  const std::vector<std::unique_ptr<ArtifactEvaluator>> evaluators_;
};

// Evaluates overdraw required to render an Artifact.
class OverdrawEvaluator : public ArtifactEvaluator {
 public:
  OverdrawEvaluator(std::unique_ptr<ingest::ViewGroupLoader> view_loader,
                    const int thread_count)
      : thread_count_(thread_count), view_loader_(std::move(view_loader)) {}
  ~OverdrawEvaluator() override = default;

  // ArtifactEvaluator implementation.
  std::string Name() const override { return "Overdraw"; }
  std::vector<Result> Evaluate(
      const artifact::Artifact& artifact) const override;

 private:
  const int thread_count_;

  // Loads the original cameras.
  std::unique_ptr<ingest::ViewGroupLoader> view_loader_;
};

// Evaluates an Artifact by simulating rendering of the Artifact to match the
// original views.
class RenderEvaluator : public ArtifactEvaluator {
 public:
  RenderEvaluator(int thread_count, std::unique_ptr<RenderSim> render_sim,
                  std::unique_ptr<ingest::ViewGroupLoader> view_loader)
      : thread_count_(thread_count),
        view_loader_(std::move(view_loader)),
        render_sim_(std::move(render_sim)) {}
  ~RenderEvaluator() override = default;

  // ArtifactEvaluator implementation.
  std::string Name() const override { return "Render"; }
  std::vector<Result> Evaluate(
      const artifact::Artifact& artifact) const override;

 private:
  const int thread_count_;

  // Loads the original views.
  std::unique_ptr<ingest::ViewGroupLoader> view_loader_;

  // Simulates rendering of the Artifact as a textured mesh.
  std::unique_ptr<RenderSim> render_sim_;
};

// Measures distortion of the original (full resolution) point cloud when
// projected onto the final Artifact's mesh.
//
// Distortion is measured by the (scaled) great-circle-distance between original
// points & their "warped" points after projection onto the nearest position
// mesh.  As a result, this is resolution agnostic.
class GeometryDistortionEvaluator : public ArtifactEvaluator {
 public:
  GeometryDistortionEvaluator(
      int thread_count, std::unique_ptr<ingest::ViewGroupLoader> view_loader)
      : thread_count_(thread_count), view_loader_(std::move(view_loader)) {}
  ~GeometryDistortionEvaluator() override = default;

  // ArtifactEvaluator implementation.
  std::string Name() const override { return "GeometryDistortion"; }
  std::vector<Result> Evaluate(
      const artifact::Artifact& artifact) const override;

 private:
  const int thread_count_;

  // Loads the original views.
  std::unique_ptr<ingest::ViewGroupLoader> view_loader_;
};

}  // namespace artifact
}  // namespace seurat

#endif  // VR_SEURAT_ARTIFACT_EVALUATION_ARTIFACT_EVALUATOR_H_
