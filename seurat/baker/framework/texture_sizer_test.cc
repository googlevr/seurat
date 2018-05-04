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

#include "seurat/baker/framework/texture_sizer.h"

#include <algorithm>
#include <numeric>

#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/baker/framework/frame.h"
#include "seurat/image/atlaser.h"
#include "seurat/image/fixed_width_atlaser.h"

namespace seurat {
namespace baker {
namespace {

using image::Atlaser;
using image::FixedWidthAtlaser;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2i;
using ion::math::Vector3f;

// Returns a dense grid with z=-1 and with (x, y) in [-1.0, 1.0]^2.
//
// This is essentially one face of a cubemap.
std::vector<Frame> MakeTessellatedQuad() {
  const int kResolution = 10;
  std::vector<Frame> frames;
  for (int y = 0; y < kResolution; ++y) {
    for (int x = 0; x < kResolution; ++x) {
      std::array<Point3f, 4> quad = {{Point3f(x, y, 0.0f),          //
                                      Point3f(x + 1, y, 0.0f),      //
                                      Point3f(x + 1, y + 1, 0.0f),  //
                                      Point3f(x, y + 1, 0.0f)}};
      for (Point3f& p : quad) {
        p /= static_cast<float>(kResolution);
        p -= {0.5f, 0.5f, 0.0f};
        p *= 2.0f;
        p[2] = -1.0f;
      }
      Frame frame;
      frame.draw_order = 0;
      frame.quad = quad;
      frames.push_back(frame);
    }
  }

  return frames;
}

TEST(TextureSizerTest, AreaTextureSizer) {
  Point3f p(50.0f, 30.0f, -100.0f);
  Vector3f v1(0.3f, -8.0f, 5.0f);
  Vector3f v2(-4.5f, -2.0f, 1.0f);

  Frame frame = {0, {{p, p + v1, p + v1 + v2, p + v2}}};

  Vector2i size;

  AreaTextureSizer sizer;
  sizer.ComputeTextureSizes(absl::Span<const Frame>(&frame, 1),
                            absl::Span<Vector2i>(&size, 1));
  EXPECT_EQ(static_cast<int>(ion::math::Length(v1)), size[0]);
  EXPECT_EQ(static_cast<int>(ion::math::Length(v2)), size[1]);
}

TEST(TextureSizerTest, ProjectedAreaTextureSizer) {
  // Position the frame very far away so that the small-angle approximation
  // holds.
  Point3f p(50.0f, 30.0f, -1000.0f);
  Vector3f v1(0.3f, -8.0f, 5.0f);
  Vector3f v2(-4.5f, -2.0f, 1.0f);

  Frame frame0 = {0, {{p, p + v1, p + v1 + v2, p + v2}}};
  // 4x the area of frame0.
  Frame frame1 = {
      0, {{p, p + v1 * 2.0f, p + v1 * 2.0f + v2 * 2.0f, p + v2 * 2.0f}}};
  std::array<Frame, 2> frames = {{frame0, frame1}};

  ProjectedAreaTextureSizer sizer(11);
  std::array<Vector2i, 2> sizes;
  sizer.ComputeTextureSizes(frames, absl::MakeSpan(sizes));
  EXPECT_LT(0, sizes[0][0]);
  EXPECT_LT(0, sizes[0][1]);
  EXPECT_LT(0, sizes[1][0]);
  EXPECT_LT(0, sizes[1][1]);

  // frame1 should be nearly double the size of frame0.
  EXPECT_NEAR(2.0f, sizes[1][0] / static_cast<float>(sizes[0][0]), 1e-2f);
  EXPECT_NEAR(2.0f, sizes[1][1] / static_cast<float>(sizes[0][1]), 1e-2f);
}

TEST(TextureSizerTest, ProjectedAreaTextureSizer_Grid) {
  // Test with a grid of known total area.
  const float kPixelsPerDegree = 10;

  std::vector<Frame> frames = MakeTessellatedQuad();
  std::vector<Vector2i> sizes(frames.size());

  ProjectedAreaTextureSizer sizer(kPixelsPerDegree);
  sizer.ComputeTextureSizes(frames, absl::MakeSpan(sizes));

  float radius = (kPixelsPerDegree * 360.0f) / (2.0f * M_PI);
  float total_resolution_of_one_sphere = 4.0f * M_PI * radius * radius;

  // The grid covers one face of a cube.
  float expected_size = total_resolution_of_one_sphere / 6.0f;

  int total_size = std::accumulate(
      sizes.begin(), sizes.end(), 0,
      [](int area, const Vector2i& size) { return area + size[0] * size[1]; });
  // Allow some slack for the small angle approximation.
  EXPECT_NEAR(1.0f, total_size / expected_size, 1.0e-2f);
}

TEST(TextureSizerTest, BucketTextureSizer) {
  Vector2i block_size(8, 5);
  BucketTextureSizer bucket_sizer(
      std::unique_ptr<TextureSizer>(new AreaTextureSizer), block_size);

  std::array<Frame, 5> frames = {
      {// A 1x1 frame.
       {0, {{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}}},
       // A 2x2 frame.
       {0, {{{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}}}},
       // A 3x2 frame.
       {0, {{{0, 0, 0}, {3, 0, 0}, {3, 2, 0}, {0, 2, 0}}}},
       // A 3x1 frame.
       {0, {{{0, 0, 0}, {3, 0, 0}, {3, 1, 0}, {0, 1, 0}}}},
       // A 99x99 frame.
       {0, {{{0, 0, 0}, {99, 0, 0}, {99, 99, 0}, {0, 99, 0}}}}}};
  std::array<Vector2i, 5> sizes;
  bucket_sizer.ComputeTextureSizes(frames, absl::MakeSpan(sizes));
  EXPECT_EQ(Vector2i(8, 5), sizes[0]);
  EXPECT_TRUE(Vector2i(8, 5) == sizes[1]);
  EXPECT_EQ(Vector2i(8, 5), sizes[2]);
  EXPECT_EQ(Vector2i(8, 5), sizes[3]);
  EXPECT_TRUE(Vector2i(104, 100) == sizes[4]);
}

TEST(TextureSizerTest, ConstrainedAtlasTextureSizer) {
  const Vector2i max_texture_size(50, 50);
  const Vector2i block_size(8, 5);
  const float kPixelsPerDegree = 10;
  std::shared_ptr<Atlaser> atlaser(new FixedWidthAtlaser(max_texture_size));
  ConstrainedAtlasTextureSizer constrained_atlas_sizer(
      atlaser, [=](float scale) {
        return std::unique_ptr<TextureSizer>(new BucketTextureSizer(
            std::unique_ptr<TextureSizer>(
                new ProjectedAreaTextureSizer(kPixelsPerDegree * scale)),
            block_size));
      });

  std::array<Frame, 5> frames = {
      {// A 1x1 frame.
       {0, {{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}}},
       // A 2x2 frame.
       {0, {{{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}}}},
       // A 3x2 frame.
       {0, {{{0, 0, 0}, {3, 0, 0}, {3, 2, 0}, {0, 2, 0}}}},
       // A 3x1 frame.
       {0, {{{0, 0, 0}, {3, 0, 0}, {3, 1, 0}, {0, 1, 0}}}},
       // A 99x99 frame.
       {0, {{{0, 0, 0}, {99, 0, 0}, {99, 99, 0}, {0, 99, 0}}}}}};

  std::array<Vector2i, 5> sizes;
  constrained_atlas_sizer.ComputeTextureSizes(frames, absl::MakeSpan(sizes));

  Vector2i total_texture_size;
  std::vector<Point2i> layout(sizes.size());
  atlaser->LayoutTiles(sizes, &total_texture_size, absl::MakeSpan(layout));
  EXPECT_GE(max_texture_size[0], total_texture_size[0]);
  EXPECT_GE(max_texture_size[1], total_texture_size[1]);
}

TEST(TextureSizerTest, ConstrainedAtlasTextureSizer_Unconstrained) {
  const Vector2i max_texture_size(std::numeric_limits<int>::max(),
                                  std::numeric_limits<int>::max());
  const Vector2i block_size(8, 5);
  const float kPixelsPerDegree = 10;
  std::shared_ptr<Atlaser> atlaser(new FixedWidthAtlaser(max_texture_size));
  ConstrainedAtlasTextureSizer constrained_atlas_sizer(
      atlaser, [=](float scale) {
        return std::unique_ptr<TextureSizer>(new BucketTextureSizer(
            std::unique_ptr<TextureSizer>(
                new ProjectedAreaTextureSizer(kPixelsPerDegree * scale)),
            block_size));
      });

  std::array<Frame, 5> frames = {
      {// A 1x1 frame.
       {0, {{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}}},
       // A 2x2 frame.
       {0, {{{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}}}},
       // A 3x2 frame.
       {0, {{{0, 0, 0}, {3, 0, 0}, {3, 2, 0}, {0, 2, 0}}}},
       // A 3x1 frame.
       {0, {{{0, 0, 0}, {3, 0, 0}, {3, 1, 0}, {0, 1, 0}}}},
       // A 99x99 frame.
       {0, {{{0, 0, 0}, {99, 0, 0}, {99, 99, 0}, {0, 99, 0}}}}}};

  std::array<Vector2i, 5> sizes;
  constrained_atlas_sizer.ComputeTextureSizes(frames, absl::MakeSpan(sizes));

  Vector2i total_texture_size;
  std::vector<Point2i> layout(sizes.size());
  atlaser->LayoutTiles(sizes, &total_texture_size, absl::MakeSpan(layout));
  EXPECT_GE(max_texture_size[0], total_texture_size[0]);
  EXPECT_GE(max_texture_size[1], total_texture_size[1]);

  std::array<Vector2i, 5> sizes_with_unit_scale;
  BucketTextureSizer(std::unique_ptr<TextureSizer>(
                         new ProjectedAreaTextureSizer(kPixelsPerDegree)),
                     block_size)
      .ComputeTextureSizes(frames, absl::MakeSpan(sizes_with_unit_scale));
  EXPECT_EQ(sizes_with_unit_scale, sizes);
}

}  // namespace
}  // namespace baker
}  // namespace seurat
