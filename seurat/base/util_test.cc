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

#include "seurat/base/util.h"

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/color.h"

namespace seurat {
namespace base {

TEST(VectorUtilsTest, ConvertVector) {
  ion::math::Vector3f ion_vec3f(4.6f, 2.8f, -8.6f);
  Color3f color_vec3f(4.6f, 2.8f, -8.6f);
  EXPECT_EQ(color_vec3f, ConvertVector<Color3f>(ion_vec3f));
  EXPECT_EQ(ion_vec3f, ConvertVector<ion::math::Vector3f>(color_vec3f));

  ion::math::Vector2i ion_vec2i(1, 2);
  Color2i color_vec2i(1, 2);
  EXPECT_EQ(color_vec2i, ConvertVector<Color2i>(ion_vec2i));
  EXPECT_EQ(ion_vec2i, ConvertVector<ion::math::Vector2i>(color_vec2i));

  ion::math::Vector4f ion_vec4f(2.0f, 4.0f, 6.0f, 8.0f);
  Color4ui16 color_vec4ui16(2, 4, 6, 8);
  EXPECT_EQ(ion_vec4f, ConvertVector<ion::math::Vector4f>(color_vec4ui16));
}

TEST(RoundUtilsTest, RoundNumbers) {
  EXPECT_EQ(RoundModN(6, 7), 7);
  EXPECT_EQ(RoundModN(7, 7), 7);
  EXPECT_EQ(RoundModN(8, 7), 14);
  for (int mod = 1; mod <= 10; ++mod) {
    for (int x = -10; x <= 10; ++x) {
      int x_rounded = RoundModN(x, mod);
      int x_r1 = ((x + mod - 1) / mod) * mod;
      EXPECT_EQ(x_rounded, x_r1);
    }
  }
}

}  // namespace base
}  // namespace seurat
