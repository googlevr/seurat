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

#include "seurat/geometry/binning_point_cloud_builder.h"

#include <vector>

#include "ion/math/matrixutils.h"
#include "ion/math/range.h"
#include "ion/math/transformutils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/geometry/cube_face.h"

namespace seurat {
namespace geometry {
namespace {

using geometry::CubeFace;
using ion::math::Matrix4f;
using ion::math::Point3f;
using ion::math::Vector3i;

TEST(BinningPointCloudBuilder, ComputeResolutionAndNearClip) {
  const ion::math::Range3f kHeadbox({-0.1f, -0.2f, -0.3f}, {0.1f, 0.2f, 0.3f});
  const float kPixelsPerDegree = 1.0f;
  Vector3i resolution;
  float near_clip;
  BinningPointCloudBuilder::ComputeResolutionAndNearClip(
      kHeadbox, kPixelsPerDegree, &resolution, &near_clip);
  EXPECT_EQ(resolution, Vector3i(180, 180, 90));
  EXPECT_EQ(near_clip, 0.1f);
}

TEST(BinningPointCloudBuilder, Build) {
  const int kThreadCount = 8;
  const Vector3i kResolution(4, 8, 2);
  const int kPointsPerFace = kResolution[0] * kResolution[1] * kResolution[2];
  const float kNear = 0.01f;

  std::array<Matrix4f, 6> world_from_clip;
  const Matrix4f clip_from_eye(          //
      1.0f, 0.0f, 0.0f, 0.0f,            //
      0.0f, 1.0f, 0.0f, 0.0f,            //
      0.0f, 0.0f, -1.0f, -2.0f * kNear,  //
      0.0f, 0.0f, -1.0f, 0.0f);
  for (int cube_face_index = 0; cube_face_index < 6; ++cube_face_index) {
    const CubeFace cube_face = static_cast<CubeFace>(cube_face_index);
    const Matrix4f eye_from_world = LookAtMatrixFromFace(cube_face);
    const Matrix4f clip_from_world = clip_from_eye * eye_from_world;
    world_from_clip[cube_face_index] = ion::math::Inverse(clip_from_world);
  }

  BinningPointCloudBuilder builder(kThreadCount, kResolution, kNear);
  std::vector<Point3f> expected_positions;
  expected_positions.reserve(6 * kPointsPerFace);

  for (int cube_face_index = 0; cube_face_index < 6; ++cube_face_index) {
    // Generate a set of positions for each cube face that exactly matches the
    // resolution of the point cloud builder. The builder should receive exactly
    // one point per perspective grid cell.
    std::vector<Point3f> positions;
    positions.reserve(kPointsPerFace);
    for (int x = 0; x < kResolution[0]; ++x) {
      for (int y = 0; y < kResolution[1]; ++y) {
        for (int z = 0; z < kResolution[2]; ++z) {
          const Point3f position_clip(
              (x + 0.5f) / kResolution[0] * 2.0f - 1.0f,
              (y + 0.5f) / kResolution[1] * 2.0f - 1.0f,
              (z + 0.5f) / kResolution[2] * 2.0f - 1.0f);
          const Point3f position_world = ion::math::ProjectPoint(
              world_from_clip[cube_face_index], position_clip);
          positions.push_back(position_world);
          expected_positions.push_back(position_world);
        }
      }
    }
    EXPECT_TRUE(builder.AddPoints(positions).ok());
  }

  // Expect that all points were added to the point cloud.
  std::vector<Point3f> merged_positions;
  std::vector<float> weights;
  builder.GetPositionsAndWeights(&merged_positions, &weights);
  EXPECT_THAT(expected_positions,
              ::testing::UnorderedElementsAreArray(merged_positions));
  for (float weight : weights) {
    EXPECT_EQ(1.0f, weight);
  }

  // After a call to GetPositions(), the builder should be reset.
  builder.GetPositionsAndWeights(&merged_positions, &weights);
  EXPECT_EQ(0, merged_positions.size());
  EXPECT_EQ(0, weights.size());

  // Add points twice and verify that we get the same result with weights equal
  // to two.
  EXPECT_TRUE(builder.AddPoints(expected_positions).ok());
  EXPECT_TRUE(builder.AddPoints(expected_positions).ok());
  builder.GetPositionsAndWeights(&merged_positions, &weights);
  EXPECT_THAT(expected_positions,
              ::testing::UnorderedElementsAreArray(merged_positions));
  for (float weight : weights) {
    EXPECT_EQ(2.0f, weight);
  }
}

TEST(BinningPointCloudBuilder, BuildWithInputWeights) {
  const int kThreadCount = 8;
  const Vector3i kResolution(4, 8, 2);
  const int kPointsPerFace = kResolution[0] * kResolution[1] * kResolution[2];
  const float kNear = 0.01f;

  std::array<Matrix4f, 6> world_from_clip;
  const Matrix4f clip_from_eye(          //
      1.0f, 0.0f, 0.0f, 0.0f,            //
      0.0f, 1.0f, 0.0f, 0.0f,            //
      0.0f, 0.0f, -1.0f, -2.0f * kNear,  //
      0.0f, 0.0f, -1.0f, 0.0f);
  for (int cube_face_index = 0; cube_face_index < 6; ++cube_face_index) {
    const CubeFace cube_face = static_cast<CubeFace>(cube_face_index);
    const Matrix4f eye_from_world = LookAtMatrixFromFace(cube_face);
    const Matrix4f clip_from_world = clip_from_eye * eye_from_world;
    world_from_clip[cube_face_index] = ion::math::Inverse(clip_from_world);
  }

  BinningPointCloudBuilder builder(kThreadCount, kResolution, kNear);
  std::vector<Point3f> expected_positions;
  std::vector<float> expected_weights;
  expected_positions.reserve(6 * kPointsPerFace);
  expected_weights.reserve(6 * kPointsPerFace);

  for (int cube_face_index = 0; cube_face_index < 6; ++cube_face_index) {
    // Generate a set of positions for each cube face that exactly matches the
    // resolution of the point cloud builder. The builder should receive exactly
    // one point per perspective grid cell.
    std::vector<Point3f> positions;
    std::vector<float> weights;
    float weight_scale =
        1.0f / (kResolution[0] * kResolution[1] * kResolution[2]);
    positions.reserve(kPointsPerFace);
    for (int x = 0; x < kResolution[0]; ++x) {
      for (int y = 0; y < kResolution[1]; ++y) {
        for (int z = 0; z < kResolution[2]; ++z) {
          const Point3f position_clip(
              (x + 0.5f) / kResolution[0] * 2.0f - 1.0f,
              (y + 0.5f) / kResolution[1] * 2.0f - 1.0f,
              (z + 0.5f) / kResolution[2] * 2.0f - 1.0f);
          const Point3f position_world = ion::math::ProjectPoint(
              world_from_clip[cube_face_index], position_clip);
          positions.push_back(position_world);
          expected_positions.push_back(position_world);
          float weight =
              weight_scale * (kResolution[2] * (kResolution[1] * x + y) + z);
          weights.push_back(weight);
          expected_weights.push_back(weight);
        }
      }
    }
    EXPECT_TRUE(builder.AddPointsWithWeights(positions, weights).ok());
  }

  // Expect that all points were added to the point cloud.
  std::vector<Point3f> merged_positions;
  std::vector<float> merged_weights;
  builder.GetPositionsAndWeights(&merged_positions, &merged_weights);
  EXPECT_THAT(expected_positions,
              ::testing::UnorderedElementsAreArray(merged_positions));
  EXPECT_THAT(expected_weights,
              ::testing::UnorderedElementsAreArray(merged_weights));

  // After a call to GetPositions(), the builder should be reset.
  builder.GetPositionsAndWeights(&merged_positions, &merged_weights);
  EXPECT_EQ(0, merged_positions.size());
  EXPECT_EQ(0, merged_weights.size());

  // Add points twice and verify that we get the same result with doubled
  // weights.

  EXPECT_TRUE(
      builder.AddPointsWithWeights(expected_positions, expected_weights).ok());
  EXPECT_TRUE(
      builder.AddPointsWithWeights(expected_positions, expected_weights).ok());
  std::transform(expected_weights.begin(), expected_weights.end(),
                 expected_weights.begin(),
                 [](float w) -> float { return 2 * w; });

  builder.GetPositionsAndWeights(&merged_positions, &merged_weights);
  EXPECT_THAT(expected_positions,
              ::testing::UnorderedElementsAreArray(merged_positions));
  EXPECT_THAT(expected_weights,
              ::testing::UnorderedElementsAreArray(merged_weights));
}

TEST(BinningPointCloudBuilder, BuildWithHigherResolution) {
  const int kThreadCount = 8;
  const Vector3i kResolution(4, 8, 2);
  const float kNear = 0.01f;

  std::array<Matrix4f, 6> world_from_clip;
  const Matrix4f clip_from_eye(          //
      1.0f, 0.0f, 0.0f, 0.0f,            //
      0.0f, 1.0f, 0.0f, 0.0f,            //
      0.0f, 0.0f, -1.0f, -2.0f * kNear,  //
      0.0f, 0.0f, -1.0f, 0.0f);
  for (int cube_face_index = 0; cube_face_index < 6; ++cube_face_index) {
    const CubeFace cube_face = static_cast<CubeFace>(cube_face_index);
    const Matrix4f eye_from_world = LookAtMatrixFromFace(cube_face);
    const Matrix4f clip_from_world = clip_from_eye * eye_from_world;
    world_from_clip[cube_face_index] = ion::math::Inverse(clip_from_world);
  }

  BinningPointCloudBuilder builder(kThreadCount, kResolution, kNear);

  for (int cube_face_index = 0; cube_face_index < 6; ++cube_face_index) {
    // Generate a set of positions for each cube face that has two points per
    // cell of the builder's binning grid. The builder should still output
    // exactly one representative point per perspective grid cell.
    std::vector<Point3f> positions;
    const int kPointsPerFace =
        2 * kResolution[0] * kResolution[1] * kResolution[2];
    positions.reserve(kPointsPerFace);
    for (int x = 0; x < kResolution[0]; ++x) {
      for (int y = 0; y < kResolution[1]; ++y) {
        for (int z = 0; z < kResolution[2]; ++z) {
          for (int i = 1; i <= 2; ++i) {
            const Point3f position_clip(
                (x + i / 3.0f) / kResolution[0] * 2.0f - 1.0f,
                (y + i / 3.0f) / kResolution[1] * 2.0f - 1.0f,
                (z + i / 3.0f) / kResolution[2] * 2.0f - 1.0f);
            const Point3f position_world = ion::math::ProjectPoint(
                world_from_clip[cube_face_index], position_clip);
            positions.push_back(position_world);
          }
        }
      }
    }
    EXPECT_TRUE(builder.AddPoints(positions).ok());
  }

  // Expect that all points were added to the point cloud.
  std::vector<Point3f> merged_positions;
  std::vector<float> weights;
  builder.GetPositionsAndWeights(&merged_positions, &weights);
  EXPECT_EQ(6 * kResolution[0] * kResolution[1] * kResolution[2],
            merged_positions.size());
  for (float weight : weights) {
    EXPECT_EQ(2.0f, weight);
  }
}

}  // namespace
}  // namespace geometry
}  // namespace seurat
