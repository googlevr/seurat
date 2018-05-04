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

#ifndef VR_SEURAT_BASE_CAMERA_H_
#define VR_SEURAT_BASE_CAMERA_H_

#include "ion/math/matrix.h"
#include "ion/math/vector.h"

namespace seurat {
namespace base {

// This class represents the intrinsic and extrinsic parameters of a camera. The
// camera implementation also defines a depth encoding.
class Camera {
 public:
  Camera() = default;
  virtual ~Camera() = default;

  // Returns the size of the image ("sensor") in pixels.
  virtual ion::math::Vector2i GetImageSize() const = 0;

  // Returns the world-from-eye matrix (extrinsic camera parameters).
  virtual ion::math::Matrix4f GetWorldFromEye() const = 0;

  // Returns the world-space origin of the ray through the |pixel|. |pixel| must
  // be in the range [0,width] x [0,height]. This origin is defined in an
  // implementation-specific manner.  For example, a projective camera may place
  // this at the near clipping plane.
  virtual ion::math::Point3f RayOriginFloat(
      const ion::math::Point2f& pixel) const = 0;

  // Returns the world-space ray direction of the ray through the |pixel|.
  // |pixel| must be in the range [0,width] x [0,height]. The resulting vector
  // is not normalized.
  virtual ion::math::Vector3f RayDirectionFloat(
      const ion::math::Point2f& pixel) const = 0;

  // Returns the world-space position of a point seen through the |pixel| at the
  // given |depth|. |pixel| must be in the range [0,width] x [0,height]. The
  // interpretation and valid range of |depth| is implementation-specific. For
  // example, a projective camera may define this as an eye-space Z.
  virtual ion::math::Point3f RayEndFloat(const ion::math::Point2f& pixel,
                                         float depth) const = 0;

  // Returns the ray origin for an integer pixel. |pixel| must be in the range
  // [0,width) x [0,height).
  ion::math::Point3f RayOrigin(const ion::math::Point2i& pixel) const {
    return RayOriginFloat(ion::math::Point2f(pixel) +
                          ion::math::Vector2f(0.5f, 0.5f));
  }

  // Returns the ray direction for an integer pixel. |pixel| must be in the
  // range [0,width) x [0,height).
  ion::math::Vector3f RayDirection(const ion::math::Point2i& pixel) const {
    return RayDirectionFloat(ion::math::Point2f(pixel) +
                             ion::math::Vector2f(0.5f, 0.5f));
  }

  // Returns the ray end for an integer pixel. |pixel| must be in the range
  // [0,width) x [0,height).
  ion::math::Point3f RayEnd(const ion::math::Point2i& pixel,
                            float depth) const {
    return RayEndFloat(
        ion::math::Point2f(pixel) + ion::math::Vector2f(0.5f, 0.5f), depth);
  }
};

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_CAMERA_H_
