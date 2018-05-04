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

#include "seurat/geometry/quad_mesh.h"

#include <vector>

#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/base/color.h"

namespace seurat {
namespace geometry {
namespace {

using base::Color4f;
using geometry::Quad3f;
using image::Image4f;
using ion::math::Point3f;

// Test default constructor.
TEST(QuadMeshTest, DefaultConstructor) {
  QuadMesh empty_quad_mesh;
  EXPECT_EQ(empty_quad_mesh.quads.size(), 0);
  EXPECT_EQ(empty_quad_mesh.textures.size(), 0);
}

// Test indexed constructor.
TEST(QuadMeshTest, IndexedConstructor) {
  const Quad3f quad{{{0.0f, 0.0f, 0.0f},
                     {1.0f, 0.0f, 0.0f},
                     {1.0f, 1.0f, 0.0f},
                     {0.0f, 1.0f, 0.0f}}};
  const std::array<float, 4> texcoord_w{{1.0f, 1.0f, 1.0f, 1.0f}};
  const int kFirstTexture = 0;
  std::vector<IndexedQuad> quads;
  quads.push_back(IndexedQuad(quad, texcoord_w, kFirstTexture));
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  std::vector<Image4f> textures{{{2, 2}, kRed}};
  QuadMesh single_quad_mesh{quads, textures};
  EXPECT_EQ(single_quad_mesh.quads.size(), quads.size());
  EXPECT_EQ(single_quad_mesh.textures.size(), textures.size());
}

// Test indexed constructor catches errors.
TEST(QuadMeshTest, IndexConstructorValidation) {
  const Quad3f quad{{{0.0f, 0.0f, 0.0f},
                     {1.0f, 0.0f, 0.0f},
                     {1.0f, 1.0f, 0.0f},
                     {0.0f, 1.0f, 0.0f}}};
  const std::array<float, 4> texcoord_w{{1.0f, 1.0f, 1.0f, 1.0f}};
  const int kNonExistentTexture = -1;
  std::vector<IndexedQuad> invalid_quads;
  invalid_quads.push_back(IndexedQuad(quad, texcoord_w, kNonExistentTexture));
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  std::vector<Image4f> textures{{{2, 2}, kRed}};
  EXPECT_DEATH(QuadMesh(invalid_quads, textures), "non-existent texture");
}

}  // namespace
}  // namespace geometry
}  // namespace seurat
