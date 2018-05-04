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

#include "seurat/artifact/artifact_processor.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "absl/strings/str_cat.h"
#include "seurat/artifact/exr_exporter.h"
#include "seurat/artifact/ice_exporter.h"
#include "seurat/artifact/obj_exporter.h"
#include "seurat/artifact/png_exporter.h"
#include "seurat/base/color.h"
#include "seurat/base/file_system.h"
#include "seurat/geometry/mesh.h"
#include "seurat/geometry/mesh_util.h"
#include "seurat/image/image.h"
#include "seurat/mesh/mesh_component_util.h"
#include "seurat/testing/test_flags.h"

namespace seurat {
namespace artifact {
namespace {

using base::Color4f;
using base::FileSystem;
using geometry::Mesh;
using image::Image4f;

// Makes an artifact which can be exported by everything.
Artifact MakeArtifactWithEverything() {
  Artifact everything;
  everything.component = mesh::MeshComponentUtil::CreateCube({});
  auto cube_component = mesh::MeshComponentUtil::CreateCube({});
  everything.mesh =
      std::make_shared<Mesh>(mesh::MeshComponentUtil::ToMesh(*cube_component));
  Image4f texture(3, 4, Color4f(0.25f, 1.0f, 0.8f, 1.0f));
  texture.At(1, 1) = Color4f::Zero();
  everything.texture = std::make_shared<Image4f>(std::move(texture));
  return everything;
}

std::shared_ptr<FileSystem> MakeTestingFileSystem() {
  return std::make_shared<FileSystem>(testing::GetTestTmpdir());
}

void TestExporter(const ArtifactProcessor& exporter,
                  const std::vector<std::string>& filenames) {
  std::shared_ptr<FileSystem> file_system = MakeTestingFileSystem();

  Artifact artifact = MakeArtifactWithEverything();
  // Clear any existing file from a previous test.
  for (const auto& fname : filenames) {
    ASSERT_TRUE(file_system->SetContents(fname, "").ok());
  }

  // Write to file.
  Artifact final_artifact = artifact;
  auto status = exporter.Process(&final_artifact);
  EXPECT_TRUE(status.ok()) << status.error_message();

  for (const auto& fname : filenames) {
    std::string exported_file_contents;
    EXPECT_TRUE(file_system->GetContents(fname, &exported_file_contents).ok());
    // Verify that it is non-zero size.
    EXPECT_LT(0, exported_file_contents.size()) << fname;
  }
}

TEST(ExporterTest, TestAllExporters) {
  std::shared_ptr<FileSystem> fs = MakeTestingFileSystem();
  std::string basename = "seurat_artifact_exporter_test";
  std::vector<std::string> filenames;

  std::vector<std::shared_ptr<const ArtifactProcessor>> exporters;
  exporters.emplace_back(new IceExporter(fs, basename));
  exporters.emplace_back(new PngExporter(fs, basename));
  exporters.emplace_back(new ExrExporter(fs, basename));
  exporters.emplace_back(new ObjExporter(fs, basename));

  filenames.push_back(absl::StrCat(basename, ".ice"));
  filenames.push_back(absl::StrCat(basename, ".png"));
  filenames.push_back(absl::StrCat(basename, ".exr"));
  filenames.push_back(absl::StrCat(basename, ".obj"));

  TestExporter(ArtifactProcessorGroup(exporters), filenames);
}

}  // namespace
}  // namespace artifact
}  // namespace seurat
