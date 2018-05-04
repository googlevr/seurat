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

#include "seurat/tiler/build_partition.h"

#include "gtest/gtest.h"
#include "seurat/tiler/geometry_model.h"

namespace seurat {
namespace tiler {

using ion::math::Point3f;
using ion::math::Vector3f;

TEST(BuildPartitionTest, TestEquality) {
  // Add some point-indices and associated errors.
  //
  // There are two points with maximal error: 3 and 1.
  //
  // There are two points with minimal error: 8 and 9.
  BuildPartition partition1;
  partition1.AddPoint(8, 0.1f);
  partition1.AddPoint(3, 0.6f);
  partition1.AddPoint(0, 0.5f);

  // Same as partition 1, but with shuffled calls to AddPoint().
  //
  // Note that the points with minimal and maximal are added in different
  // order this time.
  BuildPartition partition2;
  partition2.AddPoint(0, 0.5f);

  EXPECT_NE(partition1, partition2);

  partition2.AddPoint(3, 0.6f);
  partition2.AddPoint(8, 0.1f);

  EXPECT_EQ(partition1, partition2);

  partition2.GetMutableModel()->center = Point3f(42.0f, 42.0f, 3.0f);
  partition2.GetMutableModel()->normal = Vector3f::AxisX();
  EXPECT_NE(partition1, partition2);

  partition1.GetMutableModel()->center = Point3f(42.0f, 42.0f, 3.0f);
  partition1.GetMutableModel()->normal = Vector3f::AxisX();
  EXPECT_EQ(partition1, partition2);
}

TEST(BuildPartitionTest, TestDeterminism) {
  // Add some point-indices and associated errors.
  //
  // There are two points with maximal error: 3 and 1.
  //
  // There are two points with minimal error: 8 and 9.
  BuildPartition partition1;
  partition1.AddPoint(8, 0.1f);
  partition1.AddPoint(3, 0.6f);
  partition1.AddPoint(0, 0.5f);
  partition1.AddPoint(2, 0.25f);
  partition1.AddPoint(1, 0.6f);
  partition1.AddPoint(9, 0.1f);

  // Same as partition 1, but with shuffled calls to AddPoint().
  //
  // Note that the points with minimal and maximal are added in different
  // order this time.
  BuildPartition partition2;
  partition2.AddPoint(0, 0.5f);
  partition2.AddPoint(9, 0.1f);
  partition2.AddPoint(2, 0.25f);
  partition2.AddPoint(8, 0.1f);
  partition2.AddPoint(1, 0.6f);
  partition2.AddPoint(3, 0.6f);

  partition1.Canonicalize();
  partition2.Canonicalize();

  EXPECT_EQ(partition1.GetTotalError(), partition2.GetTotalError());
  EXPECT_EQ(partition1.GetWorstFitPoint(), partition2.GetWorstFitPoint());
  EXPECT_EQ(partition1.GetBestFitPoint(), partition2.GetBestFitPoint());

  for (int i = 0; i < partition1.GetSize(); ++i) {
    EXPECT_EQ(partition1.GetPointIndices()[i], partition2.GetPointIndices()[i]);
  }
}

TEST(BuildPartitionTest, TestSettersAndGetters) {
  BuildPartition partition;

  EXPECT_TRUE(partition.Empty());
  EXPECT_EQ(0, partition.GetSize());
  EXPECT_EQ(0, partition.GetTotalError());

  partition.AddPoint(42, 0.5f);
  EXPECT_EQ(42, partition.GetWorstFitPoint());
  EXPECT_EQ(42, partition.GetBestFitPoint());

  partition.AddPoint(10, 0.01f);
  EXPECT_EQ(42, partition.GetWorstFitPoint());
  EXPECT_EQ(10, partition.GetBestFitPoint());

  EXPECT_EQ(0.5f + 0.01f, partition.GetTotalError());

  std::vector<int> points = partition.GetPointIndices();
  EXPECT_EQ(2, points.size());
  EXPECT_EQ(1, std::count(points.begin(), points.end(), 10));
  EXPECT_EQ(1, std::count(points.begin(), points.end(), 42));

  partition.Clear();
  EXPECT_TRUE(partition.Empty());
  EXPECT_EQ(0, partition.GetSize());
  EXPECT_EQ(0, partition.GetTotalError());
}

}  // namespace tiler
}  // namespace seurat
