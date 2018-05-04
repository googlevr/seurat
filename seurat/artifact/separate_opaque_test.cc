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

#include "seurat/artifact/separate_opaque.h"

#include <memory>

#include "ion/gfx/image.h"
#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "absl/types/span.h"
#include "seurat/artifact/artifact.h"
#include "seurat/artifact/artifact_test_util.h"
#include "seurat/base/color.h"
#include "seurat/image/image.h"
#include "seurat/image/image_util.h"
#include "seurat/image/nearly_square_atlaser.h"

namespace seurat {
namespace artifact {
namespace {

using base::Color4f;
using geometry::QuadMesh;
using image::Image4f;
using ion::math::Vector2i;

constexpr float kAlphaThreshold = 255.0f / 256.0f;

Artifact BuildInitialArtifact() {
  const Color4f kOpaqueRed{125.5f / 255.0f, 0.0f, 0.0f, 1.0f};
  const Color4f kTranslucentRed{125.5f / 255.0f, 0.0f, 0.0f, 200.0f / 255.0f};
  auto initial_quad_mesh = std::make_shared<QuadMesh>();
  *initial_quad_mesh =
      MakeMultipleQuadMesh({kOpaqueRed, kOpaqueRed, kTranslucentRed});
  Artifact initial_artifact;
  initial_artifact.quad_mesh = initial_quad_mesh;
  return initial_artifact;
}

TEST(SeparateOpaque, PassOpaque) {
  Artifact initial_artifact = BuildInitialArtifact();
  SeparateOpaque artifact_processor(SeparateOpaque::Retain::kRetainOpaque,
                                    kAlphaThreshold);
  Artifact final_artifact = initial_artifact;
  EXPECT_TRUE(artifact_processor.Process(&final_artifact).ok());
  EXPECT_TRUE(final_artifact.quad_mesh);
  EXPECT_EQ(2, final_artifact.quad_mesh->quads.size());
  EXPECT_EQ(2, final_artifact.quad_mesh->textures.size());
}

TEST(SeparateOpaque, PassTranslucent) {
  Artifact initial_artifact = BuildInitialArtifact();
  SeparateOpaque artifact_processor(SeparateOpaque::Retain::kRetainTranslucent,
                                    kAlphaThreshold);
  Artifact final_artifact = initial_artifact;
  EXPECT_TRUE(artifact_processor.Process(&final_artifact).ok());
  EXPECT_TRUE(final_artifact.quad_mesh);
  EXPECT_EQ(1, final_artifact.quad_mesh->quads.size());
  EXPECT_EQ(1, final_artifact.quad_mesh->textures.size());
}

}  // namespace
}  // namespace artifact
}  // namespace seurat
