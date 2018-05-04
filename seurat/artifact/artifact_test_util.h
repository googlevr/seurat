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

#ifndef VR_SEURAT_ARTIFACT_ARTIFACT_TEST_UTIL_H_
#define VR_SEURAT_ARTIFACT_ARTIFACT_TEST_UTIL_H_

#include <memory>

#include "absl/types/span.h"
#include "seurat/artifact/artifact.h"
#include "seurat/artifact/artifact_processor.h"
#include "seurat/base/color.h"
#include "seurat/base/status_util.h"
#include "seurat/geometry/quad_mesh.h"
#include "seurat/image/image.h"

namespace seurat {
namespace artifact {

// Returns a simple quad mesh containing a single unit quad with a solid texture
// of |quad_color|.
geometry::QuadMesh MakeSingleQuadMesh(const base::Color4f& quad_color);

// Returns a quad mesh containing two quads with a solid texture of
// |vertex_color|.
geometry::QuadMesh MakeTwoQuadMesh(const base::Color4f& vertex_color);

// Returns a quad mesh containing |quad_color.size()| quads, with the i-th quad
// having a solid texture of color |quad_colors[i]|.
geometry::QuadMesh MakeMultipleQuadMesh(
    const absl::Span<const base::Color4f> quad_colors);

// Wraps a functor to produce a particular output artifact for test.
class FakePipeline : public ArtifactProcessor {
 public:
  explicit FakePipeline(std::function<Artifact()> build_artifact)
      : build_artifact_(std::move(build_artifact)) {}

  // ArtifactProcessor implementation
  base::Status Process(Artifact* artifact) const override {
    *artifact = build_artifact_();
    return base::OkStatus();
  }

 private:
  // Builds the pipeline's output artifact.
  std::function<Artifact()> build_artifact_;
};

}  // namespace artifact
}  // namespace seurat

#endif  // VR_SEURAT_ARTIFACT_ARTIFACT_TEST_UTIL_H_
