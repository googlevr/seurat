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

#include "seurat/artifact/artifact_test_util.h"

#include <vector>

#include "ion/math/vector.h"
#include "seurat/artifact/artifact.h"
#include "seurat/base/color.h"
#include "seurat/base/util.h"
#include "seurat/image/image.h"

namespace seurat {
namespace artifact {

using base::Color4f;
using geometry::IndexedQuad;
using geometry::QuadMesh;
using image::Image4f;
using ion::math::Point3f;

QuadMesh MakeSingleQuadMesh(const Color4f& quad_color) {
  const Point3f p0{0.0f, 0.0f, 0.0f};
  const Point3f p1{1.0f, 0.0f, 0.0f};
  const Point3f p2{1.0f, 1.0f, 0.0f};
  const Point3f p3{0.0f, 1.0f, 0.0f};
  const int kFirstTexture = 0;
  std::vector<IndexedQuad> quads{
      {{{p0, p1, p2, p3}}, {{1.0f, 1.0f, 1.0f, 1.0f}}, kFirstTexture}};
  std::vector<Image4f> textures{{{2, 2}, quad_color}};
  QuadMesh single_quad_mesh{quads, textures};
  return single_quad_mesh;
}

QuadMesh MakeTwoQuadMesh(const Color4f& quad_color) {
  const Point3f p0{0.0f, 0.0f, 0.0f};
  const Point3f p1{1.0f, 0.0f, 0.0f};
  const Point3f p2{1.0f, 1.0f, 0.0f};
  const Point3f p3{0.0f, 1.0f, 0.0f};
  const int kFirstTexture = 0;
  const int kSecondTexture = 1;
  std::vector<geometry::IndexedQuad> quads{
      {{{p0, p1, p2, p3}}, {{1.0f, 1.0f, 1.0f, 1.0f}}, kFirstTexture},
      {{{p0, p1, p2, p3}}, {{1.0f, 1.0f, 1.0f, 1.0f}}, kSecondTexture}};
  std::vector<Image4f> textures{{{2, 2}, quad_color}, {{2, 4}, quad_color}};
  return QuadMesh{quads, textures};
}

QuadMesh MakeMultipleQuadMesh(absl::Span<const Color4f> quad_colors) {
  const Point3f p0{0.0f, 0.0f, 0.0f};
  const Point3f p1{1.0f, 0.0f, 0.0f};
  const Point3f p2{1.0f, 1.0f, 0.0f};
  const Point3f p3{0.0f, 1.0f, 0.0f};
  const int n_quads = quad_colors.size();
  std::vector<geometry::IndexedQuad> quads(n_quads);
  std::vector<Image4f> textures(n_quads);
  for (int i = 0; i < n_quads; ++i) {
    quads[i] = {{{p0, p1, p2, p3}}, {{1.0f, 1.0f, 1.0f, 1.0f}}, i};
    textures[i] = {{2 * i, 4 * i}, quad_colors[i]};
  }

  return QuadMesh{quads, textures};
}

}  // namespace artifact
}  // namespace seurat
