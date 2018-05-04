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

#include "seurat/baker/framework/frame_sorter.h"

#include <vector>

#include "ion/math/matrix.h"
#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/geometry/quad.h"

namespace seurat {
namespace baker {
namespace {

using geometry::Plane3f;
using geometry::Quad3f;
using ion::math::Point3f;
using ion::math::Vector3f;

Frame MakeQuadFrame(float z_depth, float half_side_length) {
  const int kInitialDrawOrder = -42;
  Quad3f quad = {{Point3f(-half_side_length, -half_side_length, z_depth),
                  Point3f(half_side_length, -half_side_length, z_depth),
                  Point3f(half_side_length, half_side_length, z_depth),
                  Point3f(-half_side_length, half_side_length, z_depth)}};
  return Frame{kInitialDrawOrder, quad};
}

TEST(FrameSorterTest, TestComputeDrawOrder_Simple) {
  // Set up a simple scenario with 3 fronto-parallel quads and points nearest to
  // their assigned quads in depth.
  const int kThreadCount = 3;

  std::vector<Frame> frames;

  // Three 2x2 quads at z=-3, z=-6, and z=-9.
  frames.push_back(MakeQuadFrame(-3.0f, 1.0f));
  frames.push_back(MakeQuadFrame(-6.0f, 1.0f));
  frames.push_back(MakeQuadFrame(-9.0f, 1.0f));

  std::vector<Point3f> points;
  // Points to be assigned to the first frame.
  points.push_back({0.0f, 0.0f, -2.9f});
  points.push_back({0.0f, 0.0f, -3.0f});
  points.push_back({0.0f, 0.0f, -3.1f});
  // Points to be assigned to the last frame.
  points.push_back({0.0f, 0.0f, -8.0f});
  points.push_back({0.0f, 0.0f, -9.0f});
  points.push_back({0.0f, 0.0f, -10.0f});
  // Points to be assigned to the middle frame.
  points.push_back({0.0f, 0.0f, -7.0f});
  points.push_back({0.0f, 0.0f, -6.0f});
  points.push_back({0.0f, 0.0f, -5.0f});

  FrameSorter frame_sorter(kThreadCount);
  frame_sorter.ComputeDrawOrder(points, absl::MakeSpan(frames));
  // Render back-to-front.
  EXPECT_EQ(0, frames[2].draw_order);
  EXPECT_EQ(1, frames[1].draw_order);
  EXPECT_EQ(2, frames[0].draw_order);
}

TEST(FrameSorterTest, TestComputeDrawOrder_Empty) {
  // Same as the "Simple" test, but with no points assigned to the first frame.
  const int kThreadCount = 3;

  std::vector<Frame> frames;

  // Three 2x2 quads at z=-3, z=-6, and z=-9.
  frames.push_back(MakeQuadFrame(-3.0f, 1.0f));
  frames.push_back(MakeQuadFrame(-6.0f, 1.0f));
  frames.push_back(MakeQuadFrame(-9.0f, 1.0f));

  std::vector<Point3f> points;
  // Points to be assigned to the last frame.
  points.push_back({0.0f, 0.0f, -8.0f});
  points.push_back({0.0f, 0.0f, -9.0f});
  points.push_back({0.0f, 0.0f, -10.0f});
  // Points to be assigned to the middle frame.
  points.push_back({0.0f, 0.0f, -7.0f});
  points.push_back({0.0f, 0.0f, -6.0f});
  points.push_back({0.0f, 0.0f, -5.0f});

  FrameSorter frame_sorter(kThreadCount);
  frame_sorter.ComputeDrawOrder(points, absl::MakeSpan(frames));
  // Render back-to-front.
  EXPECT_EQ(0, frames[2].draw_order);
  EXPECT_EQ(1, frames[1].draw_order);
  EXPECT_EQ(2, frames[0].draw_order);
}

TEST(FrameSorterTest, TestComputeDrawOrder_Reverse) {
  // Set up a simple scenario with 2 fronto-parallel quads and points which,
  // when assigned via ray intersections, result in the frames being ordered in
  // reverse.
  //
  // Note that this is not meant to be a very realistic scenario.
  const int kThreadCount = 3;

  std::vector<Frame> frames;

  // Three 2x2 quads at z=-2, and z=-10 with different sizes.
  frames.push_back(MakeQuadFrame(-2.0f, 0.1f));
  frames.push_back(MakeQuadFrame(-10.0f, 100.0f));

  std::vector<Point3f> points;
  // Points to be assigned to the first frame.
  points.push_back({0.0f, 0.0f, -4.0f});
  points.push_back({0.0f, 0.0f, -4.0f});
  points.push_back({0.0f, 0.0f, -4.0f});

  // Points to be assigned to the second frame, even though they're closer to
  // the origin.
  points.push_back({0.5f, 0.5f, -1.0});
  points.push_back({0.5f, 0.5f, -1.0});
  points.push_back({0.5f, 0.5f, -1.0});

  FrameSorter frame_sorter(kThreadCount);
  frame_sorter.ComputeDrawOrder(points, absl::MakeSpan(frames));
  // Expect front-to-back draw order.
  EXPECT_EQ(0, frames[0].draw_order);
  EXPECT_EQ(1, frames[1].draw_order);
}

}  // namespace
}  // namespace baker
}  // namespace seurat
