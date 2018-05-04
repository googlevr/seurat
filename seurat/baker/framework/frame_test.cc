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

#include "seurat/baker/framework/frame.h"

#include <algorithm>
#include <random>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace baker {
namespace {

using geometry::Plane3f;
using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Vector3f;

TEST(FrameTest, InitializeApproximateDrawOrder) {
  std::vector<Frame> frames;
  // Add a frame consisting of a quad at z=-1.
  frames.push_back({-1,
                    {{Point3f(-1.0f, -1.0f, -1.0f),  //
                      Point3f(1.0f, -1.0f, -1.0f),   //
                      Point3f(1.0f, 1.0f, -1.0f),    //
                      Point3f(-1.0f, 1.0f, -1.0f)}}});
  // Add a frame consisting of a quad at z=-3.
  frames.push_back({-1,
                    {{Point3f(-1.0f, -1.0f, -3.0f),  //
                      Point3f(1.0f, -1.0f, -3.0f),   //
                      Point3f(1.0f, 1.0f, -3.0f),    //
                      Point3f(-1.0f, 1.0f, -3.0f)}}});
  // Add a frame consisting of a quad at z=-2.
  frames.push_back({-1,
                    {{Point3f(-1.0f, -1.0f, -2.0f),  //
                      Point3f(1.0f, -1.0f, -2.0f),   //
                      Point3f(1.0f, 1.0f, -2.0f),    //
                      Point3f(-1.0f, 1.0f, -2.0f)}}});
  InitializeApproximateDrawOrder(absl::MakeSpan(frames));
  EXPECT_EQ(2, frames[0].draw_order);
  EXPECT_EQ(0, frames[1].draw_order);
  EXPECT_EQ(1, frames[2].draw_order);
}

TEST(FrameTest, PlaneFromFrame) {
  Point3f p(50.0f, 30.0f, 0.0f);
  Vector3f v1(0.3f, -8.0f, 5.0f);
  Vector3f v2(-4.5f, -2.0f, 1.0f);

  Frame frame = {0, {{p, p + v1, p + v1 + v2, p + v2}}};

  Plane3f expected_plane(p, ion::math::Cross(v1, v2));

  Plane3f actual_plane = PlaneFromFrame(frame);

  if (ion::math::Dot(expected_plane.GetNormal(), actual_plane.GetNormal()) <
      0.0f) {
    expected_plane = expected_plane.GetReversePlane();
  }

  EXPECT_NEAR(
      1.0f,
      ion::math::Dot(expected_plane.GetNormal(), actual_plane.GetNormal()),
      1e-3f);
  EXPECT_NEAR(0.0f, expected_plane.SignedDistanceToPoint(p), 1e-3f);
}

TEST(FrameTest, DilateFrame) {
  Point3f p(50.0f, 30.0f, -1000.0f);
  Vector3f v1(-4.5f, -2.0f, 1.0f);
  Vector3f v2(0.3f, -8.0f, 5.0f);
  Plane3f plane(p, ion::math::Cross(v1, v2));

  Frame frame = {
      0, {{p, p + v1, p + v1 + v2, p + v2}}, {{1.0f, 1.0f, 1.0f, 1.0f}}};

  Frame dilated = frame;
  EXPECT_TRUE(DilateFrame(100.0f, &dilated));

  // The dilated frame should still be on the same plane.
  EXPECT_NEAR(0.0f, plane.SignedDistanceToPoint(dilated.quad[0]), 1e-3f);
  EXPECT_NEAR(0.0f, plane.SignedDistanceToPoint(dilated.quad[1]), 1e-3f);
  EXPECT_NEAR(0.0f, plane.SignedDistanceToPoint(dilated.quad[2]), 1e-3f);
  EXPECT_NEAR(0.0f, plane.SignedDistanceToPoint(dilated.quad[3]), 1e-3f);

  // All of the vertices of the original frame must have valid texture
  // coordinates.
  for (int i = 0; i < 4; ++i) {
    Point3f frame_coords;
    EXPECT_TRUE(WorldToFrame(dilated, frame.quad[i], &frame_coords));
    EXPECT_LE(0.0f, frame_coords[0]);
    EXPECT_LE(0.0f, frame_coords[1]);
    EXPECT_GE(1.0f, frame_coords[0]);
    EXPECT_GE(1.0f, frame_coords[1]);
  }
}

TEST(FrameTest, DilateFrame_Degenerate) {
  Point3f p(50.0f, 30.0f, -1000.0f);

  // A degenerate quad cannot be dilated (by the current implementation), so it
  // should signal failure by returning false.
  Frame frame = {0, {{p, p, p, p}}};
  Frame dilated = frame;
  EXPECT_FALSE(DilateFrame(100.0f, &dilated));
  EXPECT_FALSE(DilateFrame(0.0f, &dilated));
}

TEST(FrameTest, WorldToFrame) {
  Point3f p(-1.0f, -1.5f, 10.0f);
  Vector3f v1(2.0f, 0.0f, 0.0f);
  Vector3f v2(0.0f, 3.0f, 0.0f);

  Frame frame = {
      0, {{p, p + v1, p + v1 + v2, p + v2}}, {{1.0f, 1.0f, 1.0f, 1.0f}}};

  Point3f frame_coords;

  // Verify that all 4 corners of the frame have the expected
  // coordinates.
  EXPECT_TRUE(WorldToFrame(frame, frame.quad[0], &frame_coords));
  EXPECT_NEAR(0.0f, ion::math::Length(frame_coords - Point3f(0.0f, 0.0f, 1.0f)),
              1e-3f);

  EXPECT_TRUE(WorldToFrame(frame, frame.quad[1], &frame_coords));
  EXPECT_NEAR(0.0f, ion::math::Length(frame_coords - Point3f(1.0f, 0.0f, 1.0f)),
              1e-3f);

  EXPECT_TRUE(WorldToFrame(frame, frame.quad[2], &frame_coords));
  EXPECT_NEAR(0.0f, ion::math::Length(frame_coords - Point3f(1.0f, 1.0f, 1.0f)),
              1e-3f);

  EXPECT_TRUE(WorldToFrame(frame, frame.quad[3], &frame_coords));
  EXPECT_NEAR(0.0f, ion::math::Length(frame_coords - Point3f(0.0f, 1.0f, 1.0f)),
              1e-3f);

  // Test with a point in the middle of the quad.
  EXPECT_TRUE(WorldToFrame(frame, p + v1 * 0.5f + v2 * 0.5f, &frame_coords));
  EXPECT_LT(0.0f, frame_coords[0]);
  EXPECT_LT(0.0f, frame_coords[1]);
  EXPECT_GT(1.0f, frame_coords[0]);
  EXPECT_GT(1.0f, frame_coords[1]);

  // Test with points outside the frame.
  EXPECT_FALSE(WorldToFrame(frame,
                            frame.quad[0]  //
                                - (frame.quad[1] - frame.quad[0]),
                            &frame_coords));
  EXPECT_FALSE(WorldToFrame(frame,
                            frame.quad[0]  //
                                - (frame.quad[3] - frame.quad[0]),
                            &frame_coords));
  EXPECT_FALSE(WorldToFrame(frame,
                            frame.quad[0]                          //
                                - (frame.quad[3] - frame.quad[0])  //
                                - (frame.quad[1] - frame.quad[0]),
                            &frame_coords));
  EXPECT_FALSE(WorldToFrame(frame,
                            frame.quad[2]  //
                                - (frame.quad[3] - frame.quad[2]),
                            &frame_coords));
  EXPECT_FALSE(WorldToFrame(frame,
                            frame.quad[2]  //
                                - (frame.quad[1] - frame.quad[2]),
                            &frame_coords));
  EXPECT_FALSE(WorldToFrame(frame,
                            frame.quad[2]                          //
                                - (frame.quad[3] - frame.quad[2])  //
                                - (frame.quad[1] - frame.quad[2]),
                            &frame_coords));
}

TEST(FrameTest, FrameToWorld) {
  const float kEpsilon = 1e-5f;

  Point3f p(-1.0f, -1.5f, 10.0f);
  Vector3f v1(2.0f, 0.0f, 0.0f);
  Vector3f v2(0.0f, 3.0f, 0.0f);

  Frame frame = {
      0, {{p, p + v1, p + v1 + v2, p + v2}}, {{1.0f, 1.0f, 1.0f, 1.0f}}};

  Point3f expected_frame_coords;
  Point3f frame_coords;
  Point3f world_coords;

  expected_frame_coords = {0.2f, 0.8f, 1.0f};
  EXPECT_TRUE(FrameToWorld(frame, expected_frame_coords, &world_coords));
  EXPECT_VECTOR_NEAR(p + v1 * 0.2f + v2 * 0.8f, world_coords, kEpsilon);

  expected_frame_coords = {0.0f, 0.0f, 1.0f};
  EXPECT_TRUE(FrameToWorld(frame, expected_frame_coords, &world_coords));
  EXPECT_TRUE(WorldToFrame(frame, world_coords, &frame_coords));
  EXPECT_VECTOR_NEAR(expected_frame_coords, frame_coords, kEpsilon);

  expected_frame_coords = {1.0f, 0.0f, 1.0f};
  EXPECT_TRUE(FrameToWorld(frame, expected_frame_coords, &world_coords));
  EXPECT_TRUE(WorldToFrame(frame, world_coords, &frame_coords));
  EXPECT_VECTOR_NEAR(expected_frame_coords, frame_coords, kEpsilon);

  expected_frame_coords = {0.0f, 1.0f, 1.0f};
  EXPECT_TRUE(FrameToWorld(frame, expected_frame_coords, &world_coords));
  EXPECT_TRUE(WorldToFrame(frame, world_coords, &frame_coords));
  EXPECT_VECTOR_NEAR(expected_frame_coords, frame_coords, kEpsilon);

  expected_frame_coords = {1.0f, 1.0f, 1.0f};
  EXPECT_TRUE(FrameToWorld(frame, expected_frame_coords, &world_coords));
  EXPECT_TRUE(WorldToFrame(frame, world_coords, &frame_coords));
  EXPECT_VECTOR_NEAR(expected_frame_coords, frame_coords, kEpsilon);

  expected_frame_coords = {0.2f, 0.8f, 1.0f};
  EXPECT_TRUE(FrameToWorld(frame, expected_frame_coords, &world_coords));
  EXPECT_TRUE(WorldToFrame(frame, world_coords, &frame_coords));
  EXPECT_VECTOR_NEAR(expected_frame_coords, frame_coords, kEpsilon);

  expected_frame_coords = {0.8f, 0.2f, 1.0f};
  EXPECT_TRUE(FrameToWorld(frame, expected_frame_coords, &world_coords));
  EXPECT_TRUE(WorldToFrame(frame, world_coords, &frame_coords));
  EXPECT_VECTOR_NEAR(expected_frame_coords, frame_coords, kEpsilon);

  // Test with coordinates outside the frame.
  expected_frame_coords = {-0.1f, 0.0f, 1.0f};
  EXPECT_FALSE(FrameToWorld(frame, expected_frame_coords, &world_coords));

  expected_frame_coords = {0.0f, -0.1f, 1.0f};
  EXPECT_FALSE(FrameToWorld(frame, expected_frame_coords, &world_coords));

  expected_frame_coords = {1.1f, 0.0f, 1.0f};
  EXPECT_FALSE(FrameToWorld(frame, expected_frame_coords, &world_coords));

  expected_frame_coords = {0.0f, 1.1f, 1.0f};
  EXPECT_FALSE(FrameToWorld(frame, expected_frame_coords, &world_coords));
}

TEST(FrameTest, FreespaceRayToFrameSpace) {
  const float kEpsilon = 1e-5f;

  Point3f p(-1.0f, -1.5f, 10.0f);
  Vector3f v1(2.0f, 0.0f, 0.0f);
  Vector3f v2(0.0f, 3.0f, 0.0f);
  Plane3f plane(p, ion::math::Cross(v1, v2));

  Frame frame = {
      0, {{p, p + v1, p + v1 + v2, p + v2}}, {{1.0f, 1.0f, 1.0f, 1.0f}}};

  Point2f frame_space;

  // An arbitrary starting point which is not in the plane of the Frame.
  Point3f ray_start = p - plane.GetNormal() * 13.0f;
  Point3f frame_center = p + v1 * 0.5f + v2 * 0.5f;

  // Test with a ray going through the center of the frame.
  EXPECT_TRUE(FreespaceRayToFrameSpace(frame, ray_start,
                                       frame_center - ray_start, &frame_space));
  EXPECT_VECTOR_NEAR(Point2f(0.5f, 0.5f), frame_space, kEpsilon);

  // Test with a ray going through the center of the frame, with an endpoint off
  // of the frame.
  //
  // For a freespace ray, the corresponding 2D point should still be the center
  // of the frame at (0.5, 0.5).
  EXPECT_TRUE(FreespaceRayToFrameSpace(
      frame, ray_start, (frame_center - ray_start) * 3.0, &frame_space));
  EXPECT_VECTOR_NEAR(Point2f(0.5f, 0.5f), frame_space, kEpsilon);
}

TEST(FrameTest, SolidRayToFrameSpace) {
  const float kEpsilon = 1e-5f;

  Point3f p(-1.0f, -1.5f, 10.0f);
  Vector3f v1(2.0f, 0.0f, 0.0f);
  Vector3f v2(0.0f, 3.0f, 0.0f);
  Plane3f plane(p, ion::math::Cross(v1, v2));

  Frame frame = {
      0, {{p, p + v1, p + v1 + v2, p + v2}}, {{1.0f, 1.0f, 1.0f, 1.0f}}};

  Point2f frame_space;

  Point3f ray_start = p - plane.GetNormal() * 13.0f;
  Point3f frame_center = p + v1 * 0.5f + v2 * 0.5f;

  // Test with a ray going through the center of the frame.
  EXPECT_TRUE(
      SolidRayToFrameSpace(frame, ray_start, frame_center, &frame_space));
  EXPECT_VECTOR_NEAR(Point2f(0.5f, 0.5f), frame_space, kEpsilon);

  // Test with a ray with an endpoint behind the center of the frame, even
  // though it does not pass through the center of the frame.
  EXPECT_TRUE(SolidRayToFrameSpace(frame, ray_start + v1 * 3.0f,
                                   frame_center * 1.5f, &frame_space));
  EXPECT_VECTOR_NEAR(Point2f(0.5f, 0.5f), frame_space, kEpsilon);
}

}  // namespace
}  // namespace baker
}  // namespace seurat
