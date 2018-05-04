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

#include <array>
#include <memory>
#include <vector>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/util.h"
#include "seurat/geometry/convex_hull2d.h"
#include "seurat/geometry/polygon.h"
#include "seurat/geometry/quad.h"
#include "seurat/image/image.h"

using ion::math::Point2f;
using ion::math::Range2f;
using ion::math::Vector2f;
using seurat::image::Image4f;

namespace {

const float kEpsilon = 1.0e-5f;

}  // namespace

namespace seurat {
namespace geometry {

std::unique_ptr<Image4f> TextureSetup() {
  std::unique_ptr<Image4f> texture(new Image4f(100, 200));
  int x_c = (texture->Width() >> 1), y_c = (texture->Height() >> 1);
  for (int y = 0; y < texture->Height(); ++y) {
    for (int x = 0; x < texture->Width(); ++x) {
      float* pixel = texture->At(x, y).Data();
      pixel[0] = pixel[1] = pixel[2] = pixel[3] = 0.0f;
      int dx = x - x_c, dy = y - y_c;
      if (dx * dx + 2 * dy * dy < 500) pixel[1] = pixel[3] = 1.0f;
    }
  }
  return texture;
}

TEST(PolygonTest, TestExtent) {
  std::vector<Point2f> poly;
  poly.emplace_back(12.0f, 3.0f);
  poly.emplace_back(6.0f, 6.0f);
  poly.emplace_back(0.0f, 3.0f);
  poly.emplace_back(6.0f, 0.0f);
  const auto& extent =
      Extent(poly, {{{1.0f, 0.0f}, {0.0f, 1.0f}}}, {0.0f, 0.0f});
  EXPECT_NEAR(12.0f, (extent.GetSize())[0], kEpsilon);
  EXPECT_NEAR(6.0f, (extent.GetSize())[1], kEpsilon);
  EXPECT_NEAR(0.0f, (extent.GetMinPoint())[0], kEpsilon);
  EXPECT_NEAR(0.0f, (extent.GetMinPoint())[1], kEpsilon);
  EXPECT_NEAR(12.0f, (extent.GetMaxPoint())[0], kEpsilon);
  EXPECT_NEAR(6.0f, (extent.GetMaxPoint())[1], kEpsilon);
}

TEST(PolygonTest, TestExtent2) {
  std::vector<Point2f> poly;
  poly.emplace_back(12.0f, 3.0f);
  poly.emplace_back(6.0f, 6.0f);
  poly.emplace_back(0.0f, 3.0f);
  poly.emplace_back(6.0f, 0.0f);
  std::array<Vector2f, 2> dirs;
  dirs[0] = ion::math::Normalized(poly[2] - poly[1]);
  dirs[1] = ion::math::Orthogonal(dirs[0]);
  const auto& extent = Extent(poly, dirs, poly[1]);
  EXPECT_NEAR(10.73313f, (extent.GetSize())[0], kEpsilon);
  EXPECT_NEAR(5.36656f, (extent.GetSize())[1], kEpsilon);
  EXPECT_NEAR(-4.02492f, (extent.GetMinPoint())[0], kEpsilon);
  EXPECT_NEAR(0.0f, (extent.GetMinPoint())[1], kEpsilon);
  EXPECT_NEAR(6.70820f, (extent.GetMaxPoint())[0], kEpsilon);
  EXPECT_NEAR(5.36656f, (extent.GetMaxPoint())[1], kEpsilon);
}

TEST(PolygonTest, TestPolygonArea) {
  std::vector<Point2f> points;
  points.emplace_back(0.0f, 0.0f);
  points.emplace_back(1.0f, 0.0f);
  points.emplace_back(0.0f, 1.0f);
  float area = PolygonArea(points);
  EXPECT_NEAR(0.5f, area, kEpsilon);
}

TEST(PolygonTest, TestComputeCPolyOBB) {
  std::vector<Point2f> cpoly;
  cpoly.emplace_back(12.0f, 3.0f);
  cpoly.emplace_back(6.0f, 6.0f);
  cpoly.emplace_back(0.0f, 3.0f);
  cpoly.emplace_back(6.0f, 0.0f);
  float obb_area;
  const auto& obb = ComputeConvexPolygonOBB(cpoly, &obb_area);
  EXPECT_EQ(4, obb.size());
  EXPECT_NEAR(57.6f, obb_area, kEpsilon);
}

TEST(PolygonTest, TestTextureBoundary) {
  const auto texture = TextureSetup();
  int opaque_texel_count;
  const auto& tex_bnd = ComputeTextureBoundary(*texture, &opaque_texel_count);
  EXPECT_EQ(1115, opaque_texel_count);
  EXPECT_EQ(132, tex_bnd.size());
}

TEST(PolygonTest, TestComputeTextureConvexHull) {
  auto texture = TextureSetup();
  int opaque_texel_count;
  const auto& boundary = ComputeTextureBoundary(*texture, &opaque_texel_count);
  const auto& hull = ComputeConvexHull(boundary);
  int hull_triangle_count = hull.size() - 2;
  float hull_area = PolygonArea(hull);
  EXPECT_EQ(1115, opaque_texel_count);
  EXPECT_EQ(24, hull.size());
  EXPECT_EQ(22, hull_triangle_count);
  EXPECT_NEAR(1151, hull_area, kEpsilon);
}

TEST(PolygonTest, TestComputeOBBTilt) {
  std::vector<Point2f> cpoly;
  cpoly.emplace_back(12.0f, 3.0f);
  cpoly.emplace_back(6.0f, 6.0f);
  cpoly.emplace_back(0.0f, 3.0f);
  cpoly.emplace_back(6.0f, 0.0f);
  float obb_area;
  const auto& obb = ComputeConvexPolygonOBB(cpoly, &obb_area);
  int pivot_index;
  float phi = ComputeOBBTilt(obb, &pivot_index);
  EXPECT_NEAR(2.03444f, phi, kEpsilon);
  EXPECT_EQ(1, pivot_index);
}

TEST(PolygonTest, TestIsCounterClockwise) {
  std::vector<Point2f> convex_ccw;
  convex_ccw.push_back({0.0f, 0.0f});
  convex_ccw.push_back({1.0f, 0.0f});
  convex_ccw.push_back({1.0f, 1.0f});
  convex_ccw.push_back({0.0f, 1.0f});
  EXPECT_TRUE(IsConvexCounterClockwise(convex_ccw));

  std::vector<Point2f> nonconvex_ccw;
  nonconvex_ccw.push_back({0.0f, 0.0f});
  nonconvex_ccw.push_back({1.0f, 0.0f});
  nonconvex_ccw.push_back({0.1f, 0.1f});
  nonconvex_ccw.push_back({0.0f, 1.0f});
  EXPECT_FALSE(IsConvexCounterClockwise(nonconvex_ccw));

  std::vector<Point2f> convex_cw;
  convex_cw.push_back({0.0f, 0.0f});
  convex_cw.push_back({0.0f, 1.0f});
  convex_cw.push_back({1.0f, 1.0f});
  convex_cw.push_back({1.0f, 0.0f});
  EXPECT_FALSE(IsConvexCounterClockwise(convex_cw));

  std::vector<Point2f> nonconvex_cw;
  nonconvex_cw.push_back({0.0f, 0.0f});
  nonconvex_cw.push_back({0.0f, 1.0f});
  nonconvex_cw.push_back({0.1f, 0.1f});
  nonconvex_cw.push_back({1.0f, 0.0f});
  EXPECT_FALSE(IsConvexCounterClockwise(nonconvex_cw));
}

}  // namespace geometry
}  // namespace seurat
