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

#include "seurat/geometry/bilinear_interpolator.h"

#include <array>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/color.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace geometry {
namespace {

using base::Color4f;
using ion::math::Point2f;

TEST(BilinearInterpolatorTest, Test_Color) {
  Color4f c0(1.0f, 0.0f, 0.0f, 0.0f);
  Color4f c1(0.0f, 1.0f, 0.0f, 0.0f);
  Color4f c2(0.0f, 0.0f, 1.0f, 0.0f);
  Color4f c3(0.0f, 0.0f, 0.0f, 1.0f);

  BilinearInterpolator<Color4f> interp(
      std::array<Color4f, 4>{{c0, c1, c2, c3}});

  EXPECT_VECTOR_NEAR(c0, interp.At(0.0f, 0.0f), 1e-4f);
  EXPECT_VECTOR_NEAR(c1, interp.At(1.0f, 0.0f), 1e-4f);
  EXPECT_VECTOR_NEAR(c2, interp.At(1.0f, 1.0f), 1e-4f);
  EXPECT_VECTOR_NEAR(c3, interp.At(0.0f, 1.0f), 1e-4f);

  EXPECT_VECTOR_NEAR(0.25f * (c0 + c1 + c2 + c3), interp.At(0.5f, 0.5f), 1e-4f);
}

}  // namespace
}  // namespace geometry
}  // namespace seurat
