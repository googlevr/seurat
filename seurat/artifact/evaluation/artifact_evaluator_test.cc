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

#include "seurat/artifact/evaluation/artifact_evaluator.h"

#include "ion/math/matrix.h"
#include "ion/math/matrixutils.h"
#include "ion/math/range.h"
#include "ion/math/rangeutils.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/artifact/artifact.h"
#include "seurat/base/file_system.h"
#include "seurat/geometry/mesh.h"
#include "seurat/geometry/mesh_util.h"
#include "seurat/image/image.h"
#include "seurat/ingest/view_group_loader_test_utils.h"
#include "seurat/mesh/mesh_component_util.h"
#include "seurat/testing/test_flags.h"

namespace seurat {
namespace artifact {
namespace {

using artifact::Artifact;
using base::Color1f;
using base::Color3f;
using base::Color3i;
using base::Color4f;
using base::FileSystem;
using geometry::Mesh;
using image::Image4f;
using ingest::FakeViewGroupLoader;
using ingest::ViewGroupLoader;
using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Range3f;
using ion::math::Vector2i;
using ion::math::Vector3f;
using ion::math::Vector3i;

// Returns a mesh consisting of an origin-centered axis-aligned box with the
// side-lengths specified by |diameter|.
Mesh GetCubeMesh(const Vector3f& diameter) {
  Mesh cube =
      mesh::MeshComponentUtil::ToMesh(*mesh::MeshComponentUtil::CreateCube({}));
  cube.TransformPositions(ion::math::ScaleMatrixH(diameter));
  return cube;
}

// Returns a textured cube.
Artifact MakeTexturedCubeArtifact() {
  Artifact artifact;

  artifact.mesh = std::make_shared<Mesh>(GetCubeMesh({8.0f, 4.0f, 16.0f}));

  Image4f texture(3, 4, Color4f(0.25f, 1.0f, 0.8f, 1.0f));
  texture.At(1, 1) = Color4f::Zero();
  artifact.texture = std::make_shared<Image4f>(std::move(texture));

  return artifact;
}

std::unique_ptr<ViewGroupLoader> MakeFakeViewLoader() {
  const int kNumViewGroups = 1;
  const Vector2i kImageSize(16, 16);
  const std::array<Color3f, 6> kFaceColors = {{
      Color3f(0.0f, 0.0f, -1.0f),
      Color3f(0.0f, 0.0f, 1.0f),
      Color3f(-1.0f, 0.0f, 0.0f),
      Color3f(1.0f, 0.0f, 0.0f),
      Color3f(0.0f, -1.0f, 0.0f),
      Color3f(0.0f, 1.0f, 0.0f),
  }};
  std::array<float, 6> kFaceDepths = {{0.1f, 0.5f, 0.3f, 0.6f, 0.4f, 0.2f}};
  std::unique_ptr<ViewGroupLoader> view_group_loader(new FakeViewGroupLoader(
      kNumViewGroups, kImageSize, kFaceColors, kFaceDepths));
  return view_group_loader;
}

// Builds an ArtifactProcessor which invokes all evaluators to test.
//
// The ArtifactProcessor writes results to the specified tmp_file.
std::unique_ptr<ArtifactProcessor> BuildExporterToTest(
    const std::shared_ptr<FileSystem>& file_system,
    const std::string& tmp_stem) {
  const int kThreadCount = 3;
  const bool kEnableZBuffer = false;
  // Construct a list of all possible evaluators.
  std::vector<std::unique_ptr<ArtifactEvaluator>> evaluators;
  evaluators.emplace_back(
      new OverdrawEvaluator(MakeFakeViewLoader(), kThreadCount));
  evaluators.emplace_back(new RenderEvaluator(
      kThreadCount,
      std::unique_ptr<RenderSim>(new RenderSim(kThreadCount, kEnableZBuffer)),
      MakeFakeViewLoader()));
  evaluators.emplace_back(
      new GeometryDistortionEvaluator(kThreadCount, MakeFakeViewLoader()));
  return std::unique_ptr<ArtifactProcessor>(
      new EvaluationExporter(file_system, tmp_stem, std::move(evaluators)));
}

TEST(ArtifactEvaluatorTest, TestRunAllEvaluators) {
  // Executes all evaluators on a textured cube to sanity check that nothing
  // breaks.
  //
  // Note that the actual output is not verified, but ensuring it does not crash
  // is the least we can do here.

  std::shared_ptr<FileSystem> file_system(
      new FileSystem(testing::GetTestTmpdir()));
  std::string tmp_stem("test_artifact_evaluator_result");

  std::unique_ptr<ArtifactProcessor> exporter =
      BuildExporterToTest(file_system, tmp_stem);

  Artifact artifact = MakeTexturedCubeArtifact();
  EXPECT_TRUE(exporter->Process(&artifact).ok());

  std::string contents;
  EXPECT_TRUE(file_system->GetContents(tmp_stem + ".eval", &contents).ok());
  EXPECT_LT(0, contents.size());
}

TEST(ArtifactEvaluatorTest, TestRunAllEvaluators_EmptyArtifact) {
  // Executes all evaluators on an Artifact which doesn't support conversion to
  // anything.
  //
  // Verify that nothing breaks, even though there's nothing to evaluate.

  std::shared_ptr<FileSystem> file_system(
      new FileSystem(testing::GetTestTmpdir()));
  std::string tmp_stem("test_artifact_evaluator_result");

  std::unique_ptr<ArtifactProcessor> exporter =
      BuildExporterToTest(file_system, tmp_stem);

  Artifact empty_artifact;
  EXPECT_TRUE(exporter->Process(&empty_artifact).ok());

  std::string contents;
  EXPECT_TRUE(file_system->GetContents(tmp_stem + ".eval", &contents).ok());
  EXPECT_LT(0, contents.size());
}

}  // namespace
}  // namespace artifact
}  // namespace seurat
