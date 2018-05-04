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

#ifndef VR_SEURAT_BASE_COLOR_H_
#define VR_SEURAT_BASE_COLOR_H_

#include "ion/math/utils.h"
#include "ion/math/vector.h"

namespace seurat {
namespace base {

// Color class.
template <int Dimension, typename T>
class Color : public ion::math::VectorBase<Dimension, T> {
 public:
  // The default constructor zero-initializes all elements.
  Color() : BaseType() {}

  // Dimension-specific constructors that are passed individual element values.
  explicit Color(T e0) : BaseType(e0) {}
  Color(T e0, T e1) : BaseType(e0, e1) {}
  Color(T e0, T e1, T e2) : BaseType(e0, e1, e2) {}
  Color(T e0, T e1, T e2, T e3) : BaseType(e0, e1, e2, e3) {}

  // Constructor for a Color of dimension N from a Color of dimension N-1 and
  // a scalar of the correct type, assuming N is at least 2. This is used for
  // constructing an RGBA color from an RGB color and a separate alpha value.
  Color(const Color<Dimension - 1, T>& c, T s) : BaseType(c, s) {}

  // Constructor for a Color which copies elements from the array starting at
  // |es|.
  explicit Color(const T* es) {
    for (int c = 0; c < Dimension; ++c) {
      (*this)[c] = es[c];
    }
  }

  // Copy constructor from a Color of the same Dimension and any value type
  // that is compatible (via static_cast) with this Color's type.
  template <typename U>
  explicit Color(const Color<Dimension, U>& c) : BaseType(c) {}

  // Returns a Color containing all zeroes.
  static const Color Zero() { return ToColor(BaseType::Zero()); }

  // Returns a Color with all elements set to the given value.
  static const Color Fill(T value) { return ToColor(BaseType::Fill(value)); }

  // Convert the encoding - as opposed to solely casting the data type - from a
  // normalized float encoding to 8-bit fixed-point (i.e. 0.8) encoding.
  Color<Dimension, uint8> AsColorUI8() const {
    Color<Dimension, uint8> result;
    for (int d = 0; d < Dimension; ++d) {
      result[d] = FtoU8(BaseType::operator[](d));
    }
    return result;
  }

  // Convert the vector's data type and encoding - as opposed to solely casting
  // the data type - from a normalized 8-bit fixed-point (i.e 0.8) encoding.
  Color<Dimension, float> AsColorF() const {
    return Color<Dimension, float>(*this) / 255.0f;
  }

  // Adds another color component-wise to this.
  void operator+=(const Color& c) { BaseType::Add(c); }

  // Subtracts another color component-wise from this.
  void operator-=(const Color& c) { BaseType::Subtract(c); }

  // Component-wise mutiplies this color by another color.
  void operator*=(const Color& c) {
    *this = ToColor(BaseType::Product(*this, c));
  }

  // Component-wise divides this color by another color.
  void operator/=(const Color& c) {
    *this = ToColor(BaseType::Quotient(*this, c));
  }

  // Multiplies every channel of this color by a scalar value.
  void operator*=(T s) { BaseType::Multiply(s); }

  // Divides every channel of this color by a scalar value.
  void operator/=(T s) { BaseType::Divide(s); }

  // Unary negation.
  Color operator-() const { return ToColor(BaseType::Negation()); }

  // Returns the sum of two colors.
  friend Color operator+(const Color& c0, const Color& c1) {
    return ToColor(BaseType::Sum(c0, c1));
  }

  // Returns the difference between two colors.
  friend Color operator-(const Color& c0, const Color& c1) {
    return ToColor(BaseType::Difference(c0, c1));
  }

  // Returns a color scaled by a scalar value.
  friend Color operator*(const Color& c, T s) {
    return ToColor(BaseType::Scale(c, s));
  }

  // Returns a color scaled by a scalar value.
  friend Color operator*(T s, const Color& c) {
    // Assume the commutative property holds for type T.
    return ToColor(BaseType::Scale(c, s));
  }

  // Returns the component-wise product of two colors.
  friend Color operator*(const Color& c, const Color& s) {
    return ToColor(BaseType::Product(c, s));
  }

  // Returns a color divided by a scalar.
  friend Color operator/(const Color& c, T s) {
    return ToColor(BaseType::Divide(c, s));
  }

  // Returns the component-wise quotient of two colors.
  friend Color operator/(const Color& c, const Color& s) {
    return ToColor(BaseType::Quotient(c, s));
  }

  // Returns true if the colors are exactly equal.
  friend bool operator==(const Color& c0, const Color& c1) {
    return BaseType::AreValuesEqual(c0, c1);
  }

  // Returns false if the colors differ in any channel.
  friend bool operator!=(const Color& c0, const Color& c1) {
    return !BaseType::AreValuesEqual(c0, c1);
  }

  // Return a color with each component clamped to be <= |maximumScalar|.
  Color ClampToMax(T maximumScalar) const {
    Color result;
    for (int i = 0; i < Dimension; ++i) {
      result[i] = std::min(BaseType::operator[](i), maximumScalar);
    }
    return result;
  }

  // Return a color with each component clamped to be >= |minimumScalar|.
  Color ClampToMin(T minimumScalar) const {
    Color result;
    for (int i = 0; i < Dimension; ++i) {
      result[i] = std::max(BaseType::operator[](i), minimumScalar);
    }
    return result;
  }

  // Return a color with each component clamped to be in [0.0, 1.0].
  Color ClampToUnit() const { return ClampToMax(1.0f).ClampToMin(0.0f); }

 private:
  // Type this is derived from.
  typedef ion::math::VectorBase<Dimension, T> BaseType;

  // Converts a VectorBase of the correct type to a Color.
  static const Color ToColor(const BaseType& b) {
    // This is safe because Color is the same size as VectorBase and has no
    // virtual functions. It's ugly, but better than the alternatives.
    static_assert(
        sizeof(Color<Dimension, T>) == sizeof(BaseType),
        "Vector<> and Color<> must have the same size (and memory layout).");
    return *static_cast<const Color*>(&b);
  }

  // Convert a normalized float color component to an 8-bit range value.
  static uint8 FtoU8(float value) {
    CHECK(std::isfinite(value));
    // Clamp to unit.
    value = std::max(0.0f, value);
    value = std::min(1.0f, value);
    const uint8 kScaleToU8Range = 255;
    const float kScaledValue = value * kScaleToU8Range + 0.5f;
    return static_cast<uint8>(kScaledValue);
  }
};

// Returns true if all components of two Colors are equal within a tolerance.
template <int Dimension, typename T>
bool ColorsAlmostEqual(const base::Color<Dimension, T>& v0,
                       const base::Color<Dimension, T>& v1, T tolerance) {
  for (int i = 0; i < Dimension; ++i) {
    if (ion::math::Abs(v0[i] - v1[i]) > ion::math::Abs(tolerance)) return false;
  }
  return true;
}

// Prints a Color to a stream.
template <int Dimension, typename T>
std::ostream& operator<<(std::ostream& out, const Color<Dimension, T>& c) {
  c.Print(out, 'C');
  return out;
}

// Reads a Color from a stream.
template <int Dimension, typename T>
std::istream& operator>>(std::istream& in, Color<Dimension, T>& c) {
  c.template Read<'C'>(in);
  return in;
}

typedef Color<1, int8> Color1i8;
typedef Color<1, uint8> Color1ui8;
typedef Color<1, int16> Color1i16;
typedef Color<1, uint16> Color1ui16;
typedef Color<1, int32> Color1i;
typedef Color<1, uint32> Color1ui;
typedef Color<1, float> Color1f;
typedef Color<1, double> Color1d;
typedef Color<2, int8> Color2i8;
typedef Color<2, uint8> Color2ui8;
typedef Color<2, int16> Color2i16;
typedef Color<2, uint16> Color2ui16;
typedef Color<2, int32> Color2i;
typedef Color<2, uint32> Color2ui;
typedef Color<2, float> Color2f;
typedef Color<2, double> Color2d;
typedef Color<3, int8> Color3i8;
typedef Color<3, uint8> Color3ui8;
typedef Color<3, int16> Color3i16;
typedef Color<3, uint16> Color3ui16;
typedef Color<3, int32> Color3i;
typedef Color<3, uint32> Color3ui;
typedef Color<3, float> Color3f;
typedef Color<3, double> Color3d;
typedef Color<4, int8> Color4i8;
typedef Color<4, uint8> Color4ui8;
typedef Color<4, int16> Color4i16;
typedef Color<4, uint16> Color4ui16;
typedef Color<4, int32> Color4i;
typedef Color<4, uint32> Color4ui;
typedef Color<4, float> Color4f;
typedef Color<4, double> Color4d;

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_COLOR_H_
