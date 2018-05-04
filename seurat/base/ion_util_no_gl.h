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

#ifndef VR_SEURAT_BASE_ION_UTIL_NO_GL_H_
#define VR_SEURAT_BASE_ION_UTIL_NO_GL_H_

#include <string>

#include "ion/gfx/image.h"
#include "ion/gfx/shaderinputregistry.h"
#include "ion/gfx/texture.h"
#include "ion/gfxutils/shadermanager.h"
#include "ion/math/range.h"
#include "ion/math/vector.h"

// This file contains Ion utility functions, that don't make GL calls.

namespace seurat {
namespace base {

// Converts a Point<3, T> to a homogeneous point, represented as a Point<4, T>
// with a w-component equal to one.
template <typename T>
ion::math::Point<4, T> Point4FromPoint3(const ion::math::Point<3, T>& p) {
  return ion::math::Point<4, T>(p, 1.0f);
}

// Converts a homogeneous point represented as a Point<4, T> to a Point<3,
// T>. The w-component of |v| must not be zero.
template <typename T>
ion::math::Point<3, T> Point3FromPoint4(const ion::math::Point<4, T>& p) {
  DCHECK_NE(0.0f, p[3])
      << "p = " << p << " is a homogeneous 3D vector or point at infinity. This"
      << " method only works for homogeneous 3D points with w != 0.";
  const float inv_w = 1.0f / p[3];
  const ion::math::Point<3, T> result(p[0] * inv_w, p[1] * inv_w, p[2] * inv_w);
  DCHECK(std::isfinite(result[0]));
  DCHECK(std::isfinite(result[1]));
  DCHECK(std::isfinite(result[2]));
  return result;
}

// Projects an AABB using |transform|, returning the AABB in the new space.
template <typename T>
ion::math::Range<3, T> ProjectAABB(const ion::math::Matrix<4, T>& transform,
                                   const ion::math::Range<3, T>& aabb);

// Calculate the major axis for point |p|. NB ties are broken to the Z-axis,
// then Y-axis. As a consequence, front faces gain priority in classification.
inline int MajorAxisFromPosition(const ion::math::Point3f& p) {
  int major_axis = 0;
  if (std::fabs(p[1]) >= std::fabs(p[0])) {
    major_axis = 1;
  }
  if (std::fabs(p[2]) >= std::fabs(p[major_axis])) {
    major_axis = 2;
  }
  return major_axis;
}

// Converts a continuous (float) pixel coordinate to a discrete (int)
// pixel coordinate.
inline ion::math::Point2i DiscreteFromContinuousPixel(
    const ion::math::Point2f& p) {
  return ion::math::Point2i(static_cast<int>(std::floor(p[0])),
                            static_cast<int>(std::floor(p[1])));
}

// Creates an image with the given |format| and |dimensions|.
ion::gfx::ImagePtr CreateImage(ion::gfx::Image::Format format,
                               const ion::math::Vector2i& dimensions);

// Creates a texture from the given |image|.
ion::gfx::TexturePtr CreateTexture(ion::gfx::ImagePtr image);

// Compares two ion images for exact equality.
bool CompareImagesEqual(const ion::gfx::ImagePtr& lhs,
                        const ion::gfx::ImagePtr& rhs);

// Creates a shader program.
ion::gfx::ShaderProgramPtr CreateShaderProgram(
    const ion::gfxutils::ShaderManagerPtr& shader_manager,
    const ion::gfx::ShaderInputRegistryPtr& shader_registry,
    const string& shader_program_name, const string& vertex_shader_path,
    const string& fragment_shader_path);

// Converts a continuous (float or double) range to a discrete (int) range of
// pixels, that enclose the continuous range. In other words, it computes the
// range of discrete pixels, that are overlapped by the continuous range. The
// endpoints of the range are inclusive.
//
// Example: a range [0.3, 3.8] x [2.7, 4.2] will be converted to the pixel range
// [0, 3] x [2, 4]. This is treating pixels as squares in image space, where the
// pixel x, y covers the range [x, x+1) x [y, y+1).
template <int Dimension, typename T>
ion::math::Range<Dimension, int> EnclosingPixelRange(
    const ion::math::Range<Dimension, T>& range) {
  if (range.IsEmpty()) return ion::math::Range<Dimension, int>();

  typename ion::math::Range<Dimension, int>::Endpoint result_min;
  typename ion::math::Range<Dimension, int>::Endpoint result_max;
  for (int i = 0; i < Dimension; ++i) {
    result_min[i] = static_cast<int>(std::floor(range.GetMinPoint()[i]));
    result_max[i] = static_cast<int>(std::floor(range.GetMaxPoint()[i]));
  }
  return ion::math::Range<Dimension, int>(result_min, result_max);
}

// Returns a 4x4 perspective projection matrix with infinite far clip distance,
// otherwise the same as PerspectiveMatrixFromFrustum. The far clip epsilon may
// be zero, but when used for hardware clipping should typically be a small
// positive value that depends on the number of bits in the depth buffer, e.g.
// 2.4e-7f for 24-bit depth, or 6.1e-5f for 16-bit depth.
template <typename T>
ion::math::Matrix<4, T> PerspectiveMatrixFromInfiniteFrustum(
    T x_left, T x_right, T y_bottom, T y_top, T z_near, T z_far_epsilon) {
  const T zero = static_cast<T>(0);
  if (x_left == x_right || y_bottom == y_top || z_near <= zero) {
    return ion::math::Matrix<4, T>::Identity();
  }

  // For derivation, see for example:
  // Lengyel, E. "Projection Matrix Tricks." Game Developers Conference
  // Proceedings, 2007. http://www.terathon.com/gdc07_lengyel.pdf.
  const T X = (2 * z_near) / (x_right - x_left);
  const T Y = (2 * z_near) / (y_top - y_bottom);
  const T A = (x_right + x_left) / (x_right - x_left);
  const T B = (y_top + y_bottom) / (y_top - y_bottom);
  const T C = -1 + z_far_epsilon;
  const T D = (-2 + z_far_epsilon) * z_near;

  return ion::math::Matrix<4, T>(X, 0, A, 0, 0, Y, B, 0, 0, 0, C, D, 0, 0, -1,
                                 0);
}

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_ION_UTIL_NO_GL_H_
