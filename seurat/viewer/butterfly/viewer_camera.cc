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

#include "seurat/viewer/butterfly/viewer_camera.h"

#include "ion/math/matrixutils.h"
#include "ion/math/rotation.h"
#include "ion/math/transformutils.h"
#include "seurat/base/ion_util_no_gl.h"

namespace seurat {
namespace viewer {

using ion::math::Anglef;
using ion::math::Matrix4f;
using ion::math::Rotationf;
using ion::math::Vector3f;

namespace {

Rotationf RotationFromRollPitchYaw(const Anglef& roll, const Anglef& pitch,
                                   const Anglef& yaw) {
  Vector3f x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);
  return Rotationf::FromAxisAndAngle(z, roll) *
         (Rotationf::FromAxisAndAngle(x, pitch) *
          Rotationf::FromAxisAndAngle(y, yaw));
}

}  // namespace

ViewerCamera::ViewerCamera() {}

Matrix4f ViewerCamera::GetClipFromEye() const {
  const float kNearClip = 0.1f;
  const float kZFarEpsilon = 2.4e-7f;
  return base::PerspectiveMatrixFromInfiniteFrustum(
      -kNearClip, kNearClip, -kNearClip, kNearClip, kNearClip, kZFarEpsilon);
}

Matrix4f ViewerCamera::GetEyeFromWorld() const {
  Matrix4f rotation =
      ion::math::RotationMatrixH(RotationFromRollPitchYaw(roll_, pitch_, yaw_));
  Matrix4f translation =
      ion::math::TranslationMatrix(ion::math::Point3f::Zero() - position_);
  return rotation * translation;
}

void ViewerCamera::Translate(const ion::math::Vector3f& translation) {
  Matrix4f rotation =
      ion::math::RotationMatrixH(RotationFromRollPitchYaw(roll_, pitch_, yaw_));
  // Transform translation from eye-space to world-space.
  position_ += ion::math::Inverse(rotation) * translation;
}

void ViewerCamera::Rotate(const ion::math::Anglef& roll,
                          const ion::math::Anglef& pitch,
                          const ion::math::Anglef& yaw) {
  roll_ += roll;
  pitch_ += pitch;
  yaw_ += yaw;
}

void ViewerCamera::Reset() {
  position_ = ion::math::Point3f::Zero();
  roll_ = ion::math::Anglef();
  pitch_ = ion::math::Anglef();
  yaw_ = ion::math::Anglef();
}

}  // namespace viewer
}  // namespace seurat
