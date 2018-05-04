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

#include "seurat/base/braced_init_list.h"

#include <array>
#include <numeric>
#include <vector>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/array2d.h"

using ion::math::Point2f;

namespace seurat {
namespace base {

TEST(BracedInitListTest, ConvertBasicTypes) {
  ion::math::Vector3ui8 v3ui8(0, 0, 0);
  EXPECT_EQ("{0, 0, 0}", ToBracedInitList(v3ui8));

  ion::math::Point2i point(1, 2);
  EXPECT_EQ("{1, 2}", ToBracedInitList(point));

  std::vector<Point2f> polygon{Point2f(1.0f, 2.0f), Point2f(3.0f, 4.0f),
                               Point2f(5.0f, 6.0f), Point2f(7.0f, 8.0f)};
  EXPECT_EQ("{{1, 2}, {3, 4}, {5, 6}, {7, 8}}", ToBracedInitList(polygon));

  std::vector<int> v(9);
  std::iota(v.begin(), v.end(), 0);
  EXPECT_EQ("{0, 1, 2, 3, 4, 5, 6, 7, 8}", ToBracedInitList(v));

  std::array<int, 8> a;
  std::iota(a.begin(), a.end(), 0);
  EXPECT_EQ("{0, 1, 2, 3, 4, 5, 6, 7}", ToBracedInitList(a));
}

TEST(BracedInitListTest, ConvertArray2D) {
  Array2D<Color4f> image(2, 2);
  const Color4f kZero(0.0f, 0.0f, 0.0f, 0.0f);
  const Color4f kTest(1.0f, 4.0f, 9.0f, 25.0f);
  image.Fill(kZero);
  image.At(0, 1) = kTest;
  EXPECT_EQ("{{{0, 0, 0, 0}, {0, 0, 0, 0}},\n{{1, 4, 9, 25}, {0, 0, 0, 0}}}",
            ToBracedInitList(image));
}

}  // namespace base
}  // namespace seurat
