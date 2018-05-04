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

#ifndef VR_SEURAT_GEOMETRY_BILINEAR_INTERPOLATOR_H_
#define VR_SEURAT_GEOMETRY_BILINEAR_INTERPOLATOR_H_

#include <array>
#include "ion/math/vector.h"

namespace seurat {
namespace geometry {

// Performs bilinear interpolation of points, colors, depth values etc.
//
// ^
// | +y
//
// 3--------2
// |        |
// |        |
// |  p     |
// |        |
// 0--------1  -> +x
//
// The value at 'p' is bilinearly interpolated from the values at 0, 1, 2, 3.
//
// If p = (0,0) the result is values[0].
// If p = (1,0) the result is values[1].
template <typename T>
class BilinearInterpolator {
 public:
  explicit BilinearInterpolator(const std::array<T, 4>& values)
      : values_(values) {}

  // Returns the value interpolated at the given point, |p|, in the range
  // [0, 1]^2.
  T At(const ion::math::Point2f& p) const {
    T x_low = values_[1] * p[0] + values_[0] * (1.0f - p[0]);
    T x_high = values_[2] * p[0] + values_[3] * (1.0f - p[0]);
    return x_high * p[1] + x_low * (1.0f - p[1]);
  }

  // Returns the value interpolated at the given point, in the range [0, 1]^2.
  T At(float x, float y) const { return At({x, y}); }

 private:
  // The 4 vertices to interpolate.  See the diagram above.
  std::array<T, 4> values_;
};

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_BILINEAR_INTERPOLATOR_H_
