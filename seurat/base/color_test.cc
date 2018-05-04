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

#include "seurat/base/color.h"

#include <array>

#include "gtest/gtest.h"

//-----------------------------------------------------------------------------
// Color class tests.
//
// Note: Because Color is very similar to Vector, this tests only one
// Dimension/Scalar template type rather than all of them, as Vector does.
// -----------------------------------------------------------------------------

namespace seurat {
namespace base {

TEST(Color, ColorConstructor) {
  Color4d c4d(29.5, 30.6, 31.7, -32.8);
  EXPECT_EQ(29.5, c4d[0]);
  EXPECT_EQ(30.6, c4d[1]);
  EXPECT_EQ(31.7, c4d[2]);
  EXPECT_EQ(-32.8, c4d[3]);
}

TEST(Color, ColorCompositeConstructor) {
  Color4d c4d(Color3d(29.5, 30.6, 31.7), -32.8);
  EXPECT_EQ(29.5, c4d[0]);
  EXPECT_EQ(30.6, c4d[1]);
  EXPECT_EQ(31.7, c4d[2]);
  EXPECT_EQ(-32.8, c4d[3]);
}

TEST(Color, ArrayCopyingConstructor) {
  std::array<double, 4> color_values = {{1.0, 4.0, 3.0, 2.0}};
  Color4d c4d(color_values.data());
  EXPECT_EQ(1.0, c4d[0]);
  EXPECT_EQ(4.0, c4d[1]);
  EXPECT_EQ(3.0, c4d[2]);
  EXPECT_EQ(2.0, c4d[3]);
}

TEST(Color, ColorTypeConvertingConstructor) {
  // Integer to float.
  EXPECT_EQ(Color1f(12.0f), Color1f(Color1i(12)));
  // Integer to double.
  EXPECT_EQ(Color2d(12.0, -13.0), Color2d(Color2i(12, -13)));
  // Float to double.
  EXPECT_EQ(Color3d(12.0, -13.0, 14.5), Color3d(Color3f(12.0f, -13.0f, 14.5f)));
  // Double to integer.
  EXPECT_EQ(Color4i(12, -13, 14, 15),
            Color4i(Color4d(12.1, -13.1, 14.2, 15.0)));
}

TEST(Color, FloatFromNormalizedColorConversion) {
  Color3ui8 c3ui8(209, 30, 31);
  Color3f c3f = c3ui8.AsColorF();
  EXPECT_EQ(209 / 255.0f, c3f[0]);
  EXPECT_EQ(30 / 255.0f, c3f[1]);
  EXPECT_EQ(31 / 255.0f, c3f[2]);

  Color4ui8 c4ui8(0, 128, 191, 255);
  Color4f c4f = c4ui8.AsColorF();
  EXPECT_EQ(0.0f, c4f[0]);
  EXPECT_EQ(128 / 255.0f, c4f[1]);
  EXPECT_EQ(191 / 255.0f, c4f[2]);
  EXPECT_EQ(1.0f, c4f[3]);
}

TEST(Color, ColorEquals) {
  EXPECT_EQ(Color4d(14.1, -15.2, 16.3, -17.4),
            Color4d(14.1, -15.2, 16.3, -17.4));
  EXPECT_NE(Color4d(14.1, -15.2, 16.3, -17.4),
            Color4d(14.1, -15.2, 16.3, -17.41));
}

TEST(Color, ColorZero) {
  EXPECT_EQ(Color4d(0.0, 0.0, 0.0, 0.0), Color4d::Zero());
}

TEST(Color, ColorFill) {
  EXPECT_EQ(Color4d(1.2, 1.2, 1.2, 1.2), Color4d::Fill(1.2));
}

// Check that manipulating Colors through several encodings:
// * preserves [0.0f, 1.0f] range.
// * has expected precision.
TEST(Color, ColorEncodingExchangesZeroHalfAndOne) {
  const Color3f c3InitialZeroOneAndHalfF(0.0f, 0.5f, 1.0f);
  const Color3ui8 c3IntermediateU8 = c3InitialZeroOneAndHalfF.AsColorUI8();
  const Color3f c3FinalF = c3IntermediateU8.AsColorF();
  // The Y-axis float constant is the expected result of sending 0.5 through
  // 8-bit encoding and back to float. The constant has 9 decimal digits to
  // guarantee reconstruction.
  EXPECT_EQ(Color3f(0.0f, 0.5019608140f, 1.0f), c3FinalF);
}

TEST(Color, ColorConversionClampsOutOfRange) {
  const Color3f c3NegativeTest(-1.0f, 0.0f, 0.0f);
  const Color3ui8 c3NegativeU8 = c3NegativeTest.AsColorUI8();
  EXPECT_EQ(0, c3NegativeU8[0]);

  const Color3f c3OverflowTest(2.0f, 0.0f, 0.0f);
  const Color3ui8 c3OverflowU8 = c3OverflowTest.AsColorUI8();
  EXPECT_EQ(255, c3OverflowU8[0]);
}

TEST(Color, ColorClamp) {
  const Color3f c3InitialUnclampedF(-1.0f, 1.0f, 2.0f);
  const Color3f c3ClampedMax = c3InitialUnclampedF.ClampToMax(1.0f);
  const Color3f c3ClampedUnit = c3ClampedMax.ClampToUnit();
  const Color3f c3ExpectedClampedMaxF(-1.0f, 1.0f, 1.0f);
  const Color3f c3ExpectedClampedFinalF(0.0f, 1.0f, 1.0f);
  EXPECT_NE(c3InitialUnclampedF, c3ExpectedClampedFinalF);
  EXPECT_EQ(c3ExpectedClampedMaxF, c3ClampedMax);
  EXPECT_EQ(c3ExpectedClampedFinalF, c3ClampedUnit);
}

TEST(Color, ColorSet) {
  Color4d c4d = Color4d::Zero();
  c4d.Set(29.5, 30.6, 31.7, -32.8);
  EXPECT_EQ(Color4d(29.5, 30.6, 31.7, -32.8), c4d);
}

TEST(Color, ColorAssign) {
  Color4d c4d = Color4d::Zero();
  EXPECT_EQ(Color4d(29.5, 30.6, 31.7, -32.8),
            c4d = Color4d(29.5, 30.6, 31.7, -32.8));
  EXPECT_EQ(Color4d(29.5, 30.6, 31.7, -32.8), c4d);
}

TEST(Color, ColorMutate) {
  Color4d c4d = Color4d::Zero();
  c4d[0] = 29.5;
  c4d[1] = 30.6;
  c4d[2] = 31.7;
  c4d[3] = -32.8;
  EXPECT_EQ(Color4d(29.5, 30.6, 31.7, -32.8), c4d);
}

TEST(Color, ColorData) {
  Color4d c4d(29.5, 30.6, 31.7, -32.8);
  EXPECT_EQ(29.5, c4d.Data()[0]);
  EXPECT_EQ(30.6, c4d.Data()[1]);
  EXPECT_EQ(31.7, c4d.Data()[2]);
  EXPECT_EQ(-32.8, c4d.Data()[3]);
}

TEST(Color, ColorSelfModifyingMathOperators) {
  Color4d c0(1.0, 2.0, 3.0, 4.0);
  const Color4d c1(0.5, 1.0, 1.5, 2.0);

  c0 += Color4d(7.5, 9.5, 11.5, 13.5);
  EXPECT_EQ(Color4d(8.5, 11.5, 14.5, 17.5), c0);

  c0 -= Color4d(7.5, 9.5, 11.5, 13.5);
  EXPECT_EQ(Color4d(1.0, 2.0, 3.0, 4.0), c0);

  c0 *= 2.0;
  EXPECT_EQ(Color4d(2.0, 4.0, 6.0, 8.0), c0);

  c0 /= 4.0;
  EXPECT_EQ(Color4d(0.5, 1.0, 1.5, 2.0), c0);

  c0 += c0;
  EXPECT_EQ(Color4d(1.0, 2.0, 3.0, 4.0), c0);

  c0 *= c0;
  EXPECT_EQ(Color4d(1.0, 4.0, 9.0, 16.0), c0);

  c0 /= c1;
  EXPECT_EQ(Color4d(2.0, 4.0, 6.0, 8.0), c0);

  c0 -= c1;
  EXPECT_EQ(Color4d(1.5, 3.0, 4.5, 6.0), c0);
}

TEST(Color, ColorUnaryAndBinaryMathOperators) {
  const Color4d c0(1.5, 2.0, 6.5, -4.0);
  const Color4d c1(4.0, 5.5, 3.5, 7.0);

  // Negation.
  EXPECT_EQ(Color4d(-1.5, -2.0, -6.5, 4.0), -c0);
  EXPECT_EQ(Color4d(-4.0, -5.5, -3.5, -7.0), -c1);

  // Color * Scalar.
  EXPECT_EQ(Color4d(3.0, 4.0, 13.0, -8.0), c0 * 2.0);
  EXPECT_EQ(Color4d(3.0, 4.0, 13.0, -8.0), 2.0 * c0);

  // Color / Scalar.
  EXPECT_EQ(Color4d(3.0, 4.0, 13.0, -8.0), c0 / 0.5);

  // Color + Color.
  EXPECT_EQ(Color4d(5.5, 7.5, 10.0, 3.0), c0 + c1);
  EXPECT_EQ(Color4d(5.5, 7.5, 10.0, 3.0), c1 + c0);

  // Color - Color.
  EXPECT_EQ(Color4d(-2.5, -3.5, 3.0, -11.0), c0 - c1);

  // Color * Color.
  EXPECT_EQ(Color4d(6.0, 11.0, 22.75, -28.0), c0 * c1);
  EXPECT_EQ(Color4d(6.0, 11.0, 22.75, -28.0), c1 * c0);

  // Color / Color.
  EXPECT_EQ(Color4d(1.5 / 4.0, 2.0 / 5.5, 6.5 / 3.5, -4.0 / 7.0), c0 / c1);

  // Scaling.
  EXPECT_EQ(Color4d(6.0, 8.0, 26.0, -16.0), c0 * 4.0);
  EXPECT_EQ(Color4d(12.0, 16.5, 10.5, 21.0), 3.0 * c1);
  EXPECT_EQ(Color4d(0.75, 1.0, 3.25, -2.0), c0 / 2.0);
}

TEST(Color, ColorEqualityOperators) {
  EXPECT_TRUE(Color4d(1.5, 2.0, 6.5, -2.2) == Color4d(1.5, 2.0, 6.5, -2.2));
  EXPECT_FALSE(Color4d(1.5, 2.0, 6.5, -2.2) == Color4d(1.5, 2.0, 6.4, -2.2));
  EXPECT_FALSE(Color4d(1.5, 2.0, 6.5, -2.2) == Color4d(1.5, 2.1, 6.5, -2.2));
  EXPECT_FALSE(Color4d(1.5, 2.0, 6.5, -2.2) == Color4d(1.6, 2.0, 6.5, -2.2));
  EXPECT_FALSE(Color4d(1.5, 2.0, 6.5, -2.2) == Color4d(1.6, 2.0, 6.5, 2.2));

  EXPECT_FALSE(Color4d(1.5, 2.0, 6.5, -2.2) != Color4d(1.5, 2.0, 6.5, -2.2));
  EXPECT_TRUE(Color4d(1.5, 2.0, 6.5, -2.2) != Color4d(1.5, 2.0, 6.4, -2.2));
  EXPECT_TRUE(Color4d(1.5, 2.0, 6.5, -2.2) != Color4d(1.5, 2.1, 6.5, -2.2));
  EXPECT_TRUE(Color4d(1.5, 2.0, 6.5, -2.2) != Color4d(1.6, 2.0, 6.5, -2.2));
  EXPECT_TRUE(Color4d(1.5, 2.0, 6.5, -2.2) != Color4d(1.6, 2.0, 6.5, 2.2));
}

TEST(Color, Streaming) {
  std::ostringstream out;
  out << Color3d(4.5, 5.5, 6.5);
  EXPECT_EQ(std::string("C[4.5, 5.5, 6.5]"), out.str());

  {
    std::istringstream in("P[1.5, 2.5, 3.5]");
    Color3d c(0., 0., 0.);
    in >> c;
    EXPECT_EQ(Color3d(0., 0., 0.), c);
  }

  {
    std::istringstream in("C[1.5, 2.5, 3.5]");
    Color3d c(0., 0., 0.);
    in >> c;
    EXPECT_EQ(Color3d(1.5, 2.5, 3.5), c);
  }

  {
    std::istringstream in("C 1.5, 2.5, 3.5]");
    Color3d c(0., 0., 0.);
    in >> c;
    EXPECT_EQ(Color3d(0., 0., 0.), c);
  }

  {
    std::istringstream in("C[ 1.5, 2.5, 3.5");
    Color3d c(0., 0., 0.);
    in >> c;
    EXPECT_EQ(Color3d(0., 0., 0.), c);
  }

  {
    std::istringstream in("C[ 1.5 3.5]");
    Color2d c(0., 0.);
    in >> c;
    EXPECT_EQ(Color2d(0., 0.), c);
  }

  {
    std::istringstream in("C[ 1.5, 3.5 ]");
    Color2d c(0., 0.);
    in >> c;
    EXPECT_EQ(Color2d(1.5, 3.5), c);
  }
}

}  // namespace base
}  // namespace seurat
