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

#include "seurat/tiler/subset_geometry_solver.h"

#include <memory>
#include <numeric>
#include <vector>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/util.h"
#include "seurat/tiler/geometry_model.h"
#include "seurat/tiler/point_set.h"

namespace seurat {
namespace tiler {
namespace {

using ion::math::Point3f;
using ion::math::Vector3f;

class FakeGeometrySolver : public GeometrySolver {
 public:
  ~FakeGeometrySolver() override = default;

  // GeometrySolver implementation.
  void Init(const PointSet& point_set) override { point_set_ = point_set; }

  void InitializeModel(int point_index, GeometryModel* model) const override {
    model->center = Point3f(point_index, point_index, point_index);
    model->normal = Vector3f::AxisZ();
  }

  bool FitModel(absl::Span<const int> point_indices,
                GeometryModel* model) const override {
    const_cast<std::vector<int>&>(points_to_fit_) =
        std::vector<int>(point_indices.begin(), point_indices.end());
    return true;
  }

  float ComputeError(int point_index,
                     const GeometryModel& model) const override {
    return point_index;
  }

  // Expose the set of point indices passed to FitModel for testing.
  const std::vector<int>& GetPointsFromMostRecentFitProxy() const {
    return points_to_fit_;
  }

  // Expose the PointSet passed to Init for testing.
  PointSet GetPointSet() const { return point_set_; }

 private:
  std::vector<int> points_to_fit_;
  PointSet point_set_;
};

TEST(SubsetGeometrySolverTest, TestPassThroughMethods) {
  auto fake_solver = std::make_shared<FakeGeometrySolver>();

  const int kMaxPointCount = 3;
  SubsetGeometrySolver fast_fitting(kMaxPointCount, fake_solver);

  PointSet point_set{42, {}, {}, {}};
  fast_fitting.Init(point_set);
  EXPECT_EQ(42, fake_solver->GetPointSet().id);

  GeometryModel model;
  fast_fitting.InitializeModel(31, &model);
  EXPECT_EQ(Point3f(31, 31, 31), model.center);
  EXPECT_EQ(99.0f, fast_fitting.ComputeError(99, {}));
}

TEST(SubsetGeometrySolverTest, TestFitProxy) {
  auto fake_solver = std::make_shared<FakeGeometrySolver>();

  const int kMaxPointCount = 3;
  SubsetGeometrySolver fast_fitting(kMaxPointCount, fake_solver);

  std::vector<int> point_indices(9);
  std::iota(point_indices.begin(), point_indices.end(), 0);

  GeometryModel model;
  fast_fitting.FitModel(point_indices, &model);
  std::vector<int> subset = fake_solver->GetPointsFromMostRecentFitProxy();

  EXPECT_EQ(3, subset.size());
  int num_unique = std::unique(subset.begin(), subset.end()) - subset.begin();
  EXPECT_EQ(3, num_unique);
}

}  // namespace
}  // namespace tiler
}  // namespace seurat
