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

#ifndef VR_SEURAT_VIEWER_BUTTERFLY_VIEWER_CAMERA_H_
#define VR_SEURAT_VIEWER_BUTTERFLY_VIEWER_CAMERA_H_

#include "ion/math/angle.h"
#include "ion/math/matrix.h"
#include "ion/math/vector.h"

namespace seurat {
namespace viewer {

// A simple viewer camera model. The projection is hard-coded to use an aspect
// ratio of 1:1, a near-clip at 0.1 and a far clip plane at infinity.
class ViewerCamera {
 public:
  // Constructs a viewer camera.
  ViewerCamera();
  ~ViewerCamera() = default;

  // Returns the clip-from-eye matrix.
  ion::math::Matrix4f GetClipFromEye() const;

  // Returns the eye-from-world matrix.
  ion::math::Matrix4f GetEyeFromWorld() const;

  // Translates the camera by |translation| in eye-space.
  void Translate(const ion::math::Vector3f& translation);

  // Rotates the camera in eye-space by |roll|, |pitch| and |yaw|.
  void Rotate(const ion::math::Anglef& roll, const ion::math::Anglef& pitch,
              const ion::math::Anglef& yaw);

  // Resets the camera position to the world-space origin, looking down the
  // negative z-axis.
  void Reset();

 private:
  // World-space position of the camera.
  ion::math::Point3f position_;

  // Roll (rotation around z-axis).
  ion::math::Anglef roll_;

  // Pitch (rotation around x-axis).
  ion::math::Anglef pitch_;

  // Yaw (rotation around y-axis).
  ion::math::Anglef yaw_;
};

}  // namespace viewer
}  // namespace seurat

#endif  // VR_SEURAT_VIEWER_BUTTERFLY_VIEWER_CAMERA_H_
