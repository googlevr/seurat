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

#include "seurat/geometry/cube_face.h"

#include "ion/math/matrixutils.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "gtest/gtest.h"
#include "seurat/base/ion_util_no_gl.h"

namespace seurat {
namespace geometry {
namespace {

using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Vector3f;

const Point3f kXYZMinRight(1.0f, -1.0f, -1.0f);
const Point3f kXYZMaxRight(1.0f, 1.0f, 1.0f);
const Point3f kXYZMinLeft(-1.0f, -1.0f, 1.0f);
const Point3f kXYZMaxLeft(-1.0f, 1.0f, -1.0f);
const Point3f kXYZMinTop(-1.0f, 1.0f, -1.0f);
const Point3f kXYZMaxTop(1.0f, 1.0f, 1.0f);
const Point3f kXYZMinBottom(-1.0f, -1.0f, 1.0f);
const Point3f kXYZMaxBottom(1.0f, -1.0f, -1.0f);
const Point3f kXYZMinBack(1.0f, -1.0f, 1.0f);
const Point3f kXYZMaxBack(-1.0f, 1.0f, 1.0f);
const Point3f kXYZMinFront(-1.0f, -1.0f, -1.0f);
const Point3f kXYZMaxFront(1.0f, 1.0f, -1.0f);

TEST(CubeFaceTest, CubeFaceFromString) {
  EXPECT_EQ(CubeFace::kFront, CubeFaceFromString("front"));
  EXPECT_EQ(CubeFace::kBack, CubeFaceFromString("back"));
  EXPECT_EQ(CubeFace::kLeft, CubeFaceFromString("left"));
  EXPECT_EQ(CubeFace::kRight, CubeFaceFromString("right"));
  EXPECT_EQ(CubeFace::kTop, CubeFaceFromString("top"));
  EXPECT_EQ(CubeFace::kBottom, CubeFaceFromString("bottom"));

  EXPECT_DEATH(CubeFaceFromString("invalidstring"),
               "Invalid CubeFace name: invalidstring");
}

TEST(CubeFaceTest, IsCubeFaceString) {
  EXPECT_TRUE(IsCubeFaceString("front"));
  EXPECT_TRUE(IsCubeFaceString("back"));
  EXPECT_TRUE(IsCubeFaceString("left"));
  EXPECT_TRUE(IsCubeFaceString("right"));
  EXPECT_TRUE(IsCubeFaceString("top"));
  EXPECT_TRUE(IsCubeFaceString("bottom"));

  EXPECT_FALSE(IsCubeFaceString(""));
  EXPECT_FALSE(IsCubeFaceString("this is not a cube face"));
  EXPECT_FALSE(IsCubeFaceString("frontz"));
  EXPECT_FALSE(IsCubeFaceString("zfront"));
}

TEST(CubeFaceTest, StringFromCubeFace) {
  EXPECT_EQ("front", StringFromCubeFace(CubeFace::kFront));
  EXPECT_EQ("back", StringFromCubeFace(CubeFace::kBack));
  EXPECT_EQ("left", StringFromCubeFace(CubeFace::kLeft));
  EXPECT_EQ("right", StringFromCubeFace(CubeFace::kRight));
  EXPECT_EQ("top", StringFromCubeFace(CubeFace::kTop));
  EXPECT_EQ("bottom", StringFromCubeFace(CubeFace::kBottom));
}

void TestLookAtMatrixFromFace(CubeFace face, const Vector3f& original_direction,
                              const Vector3f& expected_direction) {
  // Verify that the specified direction is transformed to look down the
  // negative z-axis.
  Vector3f direction_in_eye_space =
      LookAtMatrixFromFace(face) * original_direction;
  float direction_dot_expected_dir =
      ion::math::Dot(direction_in_eye_space, expected_direction);
  EXPECT_EQ(1.0f, direction_dot_expected_dir);
}

TEST(CubeFaceTest, LookAtMatrixFromFace) {
  struct TestCase {
    Vector3f input_pos_cube;
    CubeFace expected_face;
  };

  //       Dotted (.) lines indicate vertical offset from X-Z plane.
  //
  //         +Y
  //   D         -Z
  //   .      |  /  B
  //   .      + /A
  //          |/ .
  // +--+--+--+--+--+-- +X
  //         /|
  //        / |
  //    C  /  |

  const TestCase kTestCases[] = {
      {Vector3f::AxisZ(), CubeFace::kBack},
      {-Vector3f::AxisZ(), CubeFace::kFront},
      {-Vector3f::AxisX(), CubeFace::kLeft},
      {Vector3f::AxisX(), CubeFace::kRight},
      {Vector3f::AxisY(), CubeFace::kTop},
      {-Vector3f::AxisY(), CubeFace::kBottom},
      {Vector3f(1.0f, 1.0f, 0.0f), CubeFace::kTop},      // A
      {Vector3f(1.0f, 0.0f, -2.0f), CubeFace::kFront},   // B
      {Vector3f(-1.0f, 0.0f, 2.0f), CubeFace::kBack},    // C
      {Vector3f(-3.0f, 2.0f, -1.0f), CubeFace::kLeft}};  // D
  for (const auto& test_case : kTestCases) {
    Point3f input_pos_cube = Point3f::Zero() + test_case.input_pos_cube;
    CubeFace calculated_face = CubeFaceFromPosition(input_pos_cube);
    EXPECT_EQ(StringFromCubeFace(calculated_face),
              StringFromCubeFace(test_case.expected_face))
        << input_pos_cube;

    // Now transform the input position from the bucketed face back to front
    // face, and check that the classification algorithm considers the point to
    // actually be on the front face.
    ion::math::Matrix4f transform_to_front =
        LookAtMatrixFromFace(calculated_face);
    Point3f point_on_front = transform_to_front * input_pos_cube;
    CubeFace expected_to_be_front_face = CubeFaceFromPosition(point_on_front);
    EXPECT_EQ(StringFromCubeFace(expected_to_be_front_face),
              StringFromCubeFace(CubeFace::kFront))
        << "Initial: " << input_pos_cube << " on "
        << StringFromCubeFace(calculated_face)
        << " On front: " << point_on_front << " Matrix: " << transform_to_front;
  }
}

TEST(CubeFaceTest, CubeFaceFromPosition) {
  TestLookAtMatrixFromFace(CubeFace::kRight, Vector3f::AxisX(),
                           -Vector3f::AxisZ());
  TestLookAtMatrixFromFace(CubeFace::kLeft, -Vector3f::AxisX(),
                           -Vector3f::AxisZ());
  TestLookAtMatrixFromFace(CubeFace::kTop, Vector3f::AxisY(),
                           -Vector3f::AxisZ());
  TestLookAtMatrixFromFace(CubeFace::kBottom, -Vector3f::AxisY(),
                           -Vector3f::AxisZ());
  TestLookAtMatrixFromFace(CubeFace::kFront, -Vector3f::AxisZ(),
                           -Vector3f::AxisZ());
  TestLookAtMatrixFromFace(CubeFace::kBack, Vector3f::AxisZ(),
                           -Vector3f::AxisZ());

  TestLookAtMatrixFromFace(CubeFace::kRight, Vector3f(1.0f, 1.0f, 1.0f),
                           Vector3f(-1.0f, 1.0f, -1.0f));
}

TEST(CubeFaceTest, XYZFromUV) {
  const Point2f kUVMin(0.0f, 0.0f);
  const Point2f kUVMax(1.0f, 1.0f);
  EXPECT_EQ(kXYZMinRight, XYZFromUV(kUVMin, CubeFace::kRight));
  EXPECT_EQ(kXYZMaxRight, XYZFromUV(kUVMax, CubeFace::kRight));
  EXPECT_EQ(kXYZMinLeft, XYZFromUV(kUVMin, CubeFace::kLeft));
  EXPECT_EQ(kXYZMaxLeft, XYZFromUV(kUVMax, CubeFace::kLeft));
  EXPECT_EQ(kXYZMinTop, XYZFromUV(kUVMin, CubeFace::kTop));
  EXPECT_EQ(kXYZMaxTop, XYZFromUV(kUVMax, CubeFace::kTop));
  EXPECT_EQ(kXYZMinBottom, XYZFromUV(kUVMin, CubeFace::kBottom));
  EXPECT_EQ(kXYZMaxBottom, XYZFromUV(kUVMax, CubeFace::kBottom));
  EXPECT_EQ(kXYZMinBack, XYZFromUV(kUVMin, CubeFace::kBack));
  EXPECT_EQ(kXYZMaxBack, XYZFromUV(kUVMax, CubeFace::kBack));
  EXPECT_EQ(kXYZMinFront, XYZFromUV(kUVMin, CubeFace::kFront));
  EXPECT_EQ(kXYZMaxFront, XYZFromUV(kUVMax, CubeFace::kFront));
}

TEST(CubeFaceTest, UVFromXYZ) {
  const Point2f kUVMin(0.0f, 0.0f);
  const Point2f kUVMax(1.0f, 1.0f);
  // Use points that are not on faces of the [-1,1]^3 cube.
  EXPECT_EQ(kUVMin, UVFromXYZ(kXYZMinRight * 42.0f, CubeFace::kRight));
  EXPECT_EQ(kUVMax, UVFromXYZ(kXYZMaxRight * 42.0f, CubeFace::kRight));
  EXPECT_EQ(kUVMin, UVFromXYZ(kXYZMinLeft * 42.0f, CubeFace::kLeft));
  EXPECT_EQ(kUVMax, UVFromXYZ(kXYZMaxLeft * 42.0f, CubeFace::kLeft));
  EXPECT_EQ(kUVMin, UVFromXYZ(kXYZMinTop * 42.0f, CubeFace::kTop));
  EXPECT_EQ(kUVMax, UVFromXYZ(kXYZMaxTop * 42.0f, CubeFace::kTop));
  EXPECT_EQ(kUVMin, UVFromXYZ(kXYZMinBottom * 42.0f, CubeFace::kBottom));
  EXPECT_EQ(kUVMax, UVFromXYZ(kXYZMaxBottom * 42.0f, CubeFace::kBottom));
  EXPECT_EQ(kUVMin, UVFromXYZ(kXYZMinBack * 42.0f, CubeFace::kBack));
  EXPECT_EQ(kUVMax, UVFromXYZ(kXYZMaxBack * 42.0f, CubeFace::kBack));
  EXPECT_EQ(kUVMin, UVFromXYZ(kXYZMinFront * 42.0f, CubeFace::kFront));
  EXPECT_EQ(kUVMax, UVFromXYZ(kXYZMaxFront * 42.0f, CubeFace::kFront));
}

TEST(CubeFaceTest, XYZFromUVW) {
  // Use a homogeneous UVW coordinate with an (arbitrary) w-value of 42.
  const Point3f kUVWMin(0.0f, 0.0f, 42.0f);
  const Point3f kUVWMax(42.0f, 42.0f, 42.0f);
  EXPECT_EQ(kXYZMinRight * 42.0f, XYZFromUVW(kUVWMin, CubeFace::kRight));
  EXPECT_EQ(kXYZMaxRight * 42.0f, XYZFromUVW(kUVWMax, CubeFace::kRight));
  EXPECT_EQ(kXYZMinLeft * 42.0f, XYZFromUVW(kUVWMin, CubeFace::kLeft));
  EXPECT_EQ(kXYZMaxLeft * 42.0f, XYZFromUVW(kUVWMax, CubeFace::kLeft));
  EXPECT_EQ(kXYZMinTop * 42.0f, XYZFromUVW(kUVWMin, CubeFace::kTop));
  EXPECT_EQ(kXYZMaxTop * 42.0f, XYZFromUVW(kUVWMax, CubeFace::kTop));
  EXPECT_EQ(kXYZMinBottom * 42.0f, XYZFromUVW(kUVWMin, CubeFace::kBottom));
  EXPECT_EQ(kXYZMaxBottom * 42.0f, XYZFromUVW(kUVWMax, CubeFace::kBottom));
  EXPECT_EQ(kXYZMinBack * 42.0f, XYZFromUVW(kUVWMin, CubeFace::kBack));
  EXPECT_EQ(kXYZMaxBack * 42.0f, XYZFromUVW(kUVWMax, CubeFace::kBack));
  EXPECT_EQ(kXYZMinFront * 42.0f, XYZFromUVW(kUVWMin, CubeFace::kFront));
  EXPECT_EQ(kXYZMaxFront * 42.0f, XYZFromUVW(kUVWMax, CubeFace::kFront));
}

TEST(CubeFaceTest, UVWFromXYZ) {
  // Use points that are not on faces of the [-1,1]^3 cube.
  const Point3f kUVWMin(0.0f, 0.0f, 42.0f);
  const Point3f kUVWMax(42.0f, 42.0f, 42.0f);
  EXPECT_EQ(kUVWMin, UVWFromXYZ(kXYZMinRight * 42.0f, CubeFace::kRight));
  EXPECT_EQ(kUVWMax, UVWFromXYZ(kXYZMaxRight * 42.0f, CubeFace::kRight));
  EXPECT_EQ(kUVWMin, UVWFromXYZ(kXYZMinLeft * 42.0f, CubeFace::kLeft));
  EXPECT_EQ(kUVWMax, UVWFromXYZ(kXYZMaxLeft * 42.0f, CubeFace::kLeft));
  EXPECT_EQ(kUVWMin, UVWFromXYZ(kXYZMinTop * 42.0f, CubeFace::kTop));
  EXPECT_EQ(kUVWMax, UVWFromXYZ(kXYZMaxTop * 42.0f, CubeFace::kTop));
  EXPECT_EQ(kUVWMin, UVWFromXYZ(kXYZMinBottom * 42.0f, CubeFace::kBottom));
  EXPECT_EQ(kUVWMax, UVWFromXYZ(kXYZMaxBottom * 42.0f, CubeFace::kBottom));
  EXPECT_EQ(kUVWMin, UVWFromXYZ(kXYZMinBack * 42.0f, CubeFace::kBack));
  EXPECT_EQ(kUVWMax, UVWFromXYZ(kXYZMaxBack * 42.0f, CubeFace::kBack));
  EXPECT_EQ(kUVWMin, UVWFromXYZ(kXYZMinFront * 42.0f, CubeFace::kFront));
  EXPECT_EQ(kUVWMax, UVWFromXYZ(kXYZMaxFront * 42.0f, CubeFace::kFront));
}

}  // namespace
}  // namespace geometry
}  // namespace seurat
