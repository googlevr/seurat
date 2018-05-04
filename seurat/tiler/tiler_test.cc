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

#include "seurat/tiler/tiler.h"

#include <cmath>
#include <random>

#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "gtest/gtest.h"
#include "seurat/base/parallel.h"
#include "seurat/base/util.h"
#include "seurat/tiler/point_set.h"
#include "seurat/tiler/tiler_test_util.h"

namespace seurat {
namespace tiler {
namespace {

using ion::math::Point3f;
using ion::math::Vector3f;

TEST(TilerTest, TestSelectionTiler) {
  const int kPointCount = 10000;
  const int kMaxTileCount = 200;
  const float kOverdrawFactor = 2.5f;

  TilerFactory::Parameters parameters;
  parameters.tile_count = kMaxTileCount;
  parameters.thread_count = 2;
  parameters.overdraw_factor = kOverdrawFactor;
  parameters.peak_overdraw_factor = kOverdrawFactor * 3.0f;
  parameters.min_subdivision_level = 1;
  parameters.max_subdivision_level = 2;
  parameters.fast = true;

  std::unique_ptr<Tiler> tiler = TilerFactory::CreateSelectionTiler(parameters);

  int point_set_id = 0;
  ExpectTilesUnitSphere(tiler.get(), point_set_id++, kPointCount);
  ExpectTilesSphericalCap(tiler.get(), point_set_id++, kPointCount,
                          Vector3f::AxisZ(), 1.0e-3f);
  ExpectTilesSphericalCap(tiler.get(), point_set_id++, kPointCount,
                          -Vector3f::AxisZ(), 1.0e-3f);
  // Test partial-scene support by testing with a spherical cap which is
  // truncated at 10-degrees off of the xz plane.
  ExpectTilesSphericalCap(tiler.get(), point_set_id++, kPointCount,
                          -Vector3f::AxisY(),
                          std::cos(2.0f * M_PI * 10.0f / 360.0f));
}

}  // namespace
}  // namespace tiler
}  // namespace seurat
