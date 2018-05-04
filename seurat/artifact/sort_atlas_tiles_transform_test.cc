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

#include "seurat/artifact/sort_atlas_tiles_transform.h"

#include <algorithm>
#include <memory>
#include <vector>

#include "gtest/gtest.h"
#include "seurat/base/array2d_util.h"
#include "seurat/base/color.h"
#include "seurat/geometry/quad_mesh.h"
#include "seurat/image/image.h"

namespace seurat {
namespace artifact {

using ion::math::Point2i;
using ion::math::Vector2i;

using base::Color4f;
using geometry::IndexedQuad;
using geometry::QuadMesh;
using image::Image4f;

namespace {

// Orders images by size lexicographically.
bool TileSizeLessThan(const Image4f& lhs, const Image4f& rhs) {
  const Vector2i& lhs_size = lhs.GetSize();
  const Vector2i& rhs_size = rhs.GetSize();
  return lhs_size[0] < rhs_size[0] ||
         (lhs_size[0] == rhs_size[0] && lhs_size[1] < rhs_size[1]);
}

}  // namespace

TEST(SortTextureTiles, TestSort) {
  // Builds a quad mesh artifact to be sorted.
  std::vector<Image4f> textures;
  textures.push_back(Image4f(Vector2i(1, 2), Color4f(0.0f, 1.0, 0.0f, 1.0f)));
  textures.push_back(Image4f(Vector2i(2, 2), Color4f(1.0f, 1.0, 0.0f, 1.0f)));
  textures.push_back(Image4f(Vector2i(2, 1), Color4f(0.0f, 1.0, 0.0f, 1.0f)));
  textures.push_back(Image4f(Vector2i(1, 1), Color4f(0.0f, 0.0, 1.0f, 1.0f)));
  std::vector<IndexedQuad> quads;
  quads.push_back(IndexedQuad({}, {}, 0));
  quads.push_back(IndexedQuad({}, {}, 1));
  quads.push_back(IndexedQuad({}, {}, 2));
  quads.push_back(IndexedQuad({}, {}, 3));
  const QuadMesh quad_mesh(std::move(quads), std::move(textures));

  Artifact quad_mesh_artifact;
  quad_mesh_artifact.quad_mesh = std::make_shared<QuadMesh>(quad_mesh);

  SortAtlasTilesTransform sort_tiles_transform;

  ASSERT_TRUE(sort_tiles_transform.Process(&quad_mesh_artifact).ok());
  ASSERT_TRUE(quad_mesh_artifact.quad_mesh);
  const QuadMesh& sorted_quad_mesh = *quad_mesh_artifact.quad_mesh;

  // The test can't continue if this fails, but failure doesn't necessarily
  // indicate the test code is buggy.
  ASSERT_EQ(quad_mesh.quads.size(), sorted_quad_mesh.quads.size());
  EXPECT_EQ(quad_mesh.textures.size(), sorted_quad_mesh.textures.size());

  EXPECT_TRUE(std::is_sorted(sorted_quad_mesh.textures.begin(),
                             sorted_quad_mesh.textures.end(),
                             TileSizeLessThan));

  // The code under test sorts the order of the textures, but leaves the quad
  // mesh defining the same data as before. Test that iterating the textures of
  // the quads of the unsorted and sorted mesh produces identical images.
  EXPECT_GT(quad_mesh.quads.size(), 0);
  for (int quad_index = 0; quad_index < quad_mesh.quads.size(); ++quad_index) {
    const Image4f& unsorted_texture =
        quad_mesh.textures[quad_mesh.quads[quad_index].texture_index];
    const Image4f& sorted_texture =
        sorted_quad_mesh
            .textures[sorted_quad_mesh.quads[quad_index].texture_index];

    // The test can't continue if this fails, but again failure doesn't
    // necessarily indicate the test code is buggy.
    ASSERT_EQ(unsorted_texture.GetSize(), sorted_texture.GetSize());
    if (unsorted_texture.GetSize() == sorted_texture.GetSize())
      base::SpatialForEachArrayEntry(
          unsorted_texture,
          [&sorted_texture](const Point2i& p, const Color4f& c) {
            EXPECT_EQ(c, sorted_texture.At(p)) << " at " << p;
          });
  }
}

}  // namespace artifact
}  // namespace seurat
