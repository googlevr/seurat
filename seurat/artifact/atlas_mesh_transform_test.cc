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

#include "seurat/artifact/atlas_mesh_transform.h"

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
using geometry::Mesh;
using geometry::QuadMesh;
using image::Image4f;
using image::NearlySquareAtlaser;
using ion::math::Vector2i;

TEST(AtlasMeshTransform, MeshAndAtlas) {
  // Prepares input data.
  const Color4f kSlightRed{125.5f / 255.0f, 0.0f, 0.0f, 1.0f};
  auto initial_quad_mesh = std::make_shared<QuadMesh>();
  *initial_quad_mesh = MakeTwoQuadMesh(kSlightRed);
  Artifact initial_artifact;
  initial_artifact.quad_mesh = initial_quad_mesh;

  // Builds the transform to be tested.
  auto atlaser = std::make_shared<NearlySquareAtlaser>();
  AtlasMeshTransform transform(atlaser);

  // Tests the transform.
  Artifact final_artifact = initial_artifact;
  EXPECT_TRUE(transform.Process(&final_artifact).ok());

  EXPECT_TRUE(final_artifact.mesh);
  EXPECT_TRUE(final_artifact.texture);
  EXPECT_TRUE(final_artifact.component);

  const Mesh& final_mesh = *final_artifact.mesh;
  const Image4f& final_atlas = *final_artifact.texture;

  // Validates final objects.
  EXPECT_EQ(1, final_mesh.GetTextureCount());
  EXPECT_EQ(8, final_mesh.GetVertexCount());
  EXPECT_EQ(4, final_mesh.GetTriangleCount());
  EXPECT_EQ(Vector2i(4, 4), final_atlas.GetSize());
}

}  // namespace
}  // namespace artifact
}  // namespace seurat
