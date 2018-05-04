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

#include "seurat/base/projective_camera.h"

#include <cmath>

#include "ion/math/matrixutils.h"
#include "ion/math/transformutils.h"

using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Point4f;
using ion::math::Vector2i;
using ion::math::Vector3f;

// If the 14th element of a projection matrix is smaller than this epsilon,
// it is considered an orthographic projection.
const float kOrthographicProjectionEpsilon = 1e-5f;

namespace {

ion::math::Matrix4f CreateViewportMatrix(const Vector2i& viewport_size) {
  const float w_2 = viewport_size[0] * 0.5f;
  const float h_2 = viewport_size[1] * 0.5f;
  return ion::math::Matrix4f(w_2, 0.0f, 0.0f, w_2,    //
                             0.0f, h_2, 0.0f, h_2,    //
                             0.0f, 0.0f, 0.5f, 0.5f,  //
                             0.0f, 0.0f, 0.0f, 1.0f);
}

}  // namespace

namespace seurat {
namespace base {

ProjectiveCamera::ProjectiveCamera(const ion::math::Vector2i& image_size,
                                   const Matrix4f& clip_from_eye,
                                   const Matrix4f& eye_from_world,
                                   const Matrix4f& world_from_eye,
                                   DepthType depth_type)
    : image_size_(image_size),
      clip_from_eye_(clip_from_eye),
      eye_from_clip_(ion::math::Inverse(clip_from_eye)),
      eye_from_world_(eye_from_world),
      world_from_eye_(world_from_eye),
      window_from_eye_(CreateViewportMatrix(image_size) * clip_from_eye),
      eye_from_window_(ion::math::Inverse(window_from_eye_)),
      window_from_world_(window_from_eye_ * eye_from_world_),
      world_from_window_(GetWorldFromEye() * eye_from_window_),
      depth_type_(depth_type) {}

Vector2i ProjectiveCamera::GetImageSize() const { return image_size_; }

Matrix4f ProjectiveCamera::GetWorldFromEye() const { return world_from_eye_; }

Point3f ProjectiveCamera::RayOriginFloat(const Point2f& pixel) const {
  // Return the point on the near-clip plane.
  return ion::math::ProjectPoint(world_from_window_, Point3f(pixel, 0.0f));
}

Vector3f ProjectiveCamera::RayDirectionFloat(const Point2f& pixel) const {
  // Return the vector from the camera position to the point on the near-clip
  // plane.
  return ion::math::ProjectPoint(world_from_window_, Point3f(pixel, 0.0f)) -
         GetWorldFromEye() * Point3f::Zero();
}

Point3f ProjectiveCamera::RayEndFloat(const Point2f& pixel, float depth) const {
  switch (depth_type_) {
    case DepthType::kWindowZ: {
      DCHECK_GE(depth, 0.0f);
      DCHECK_LE(depth, 1.0f);
      return ion::math::ProjectPoint(world_from_window_, Point3f(pixel, depth));
    }
    case DepthType::kEyeZ: {
      DCHECK_GE(depth, 0.0f);
      DCHECK(std::isfinite(depth));
      const Vector3f ray_direction_eye =
          ion::math::ProjectPoint(eye_from_window_, Point3f(pixel, 0.0f)) -
          Point3f::Zero();
      // |depth| is in the range [0, inf). Negate the value, because OpenGL
      // cameras are looking down the negative Z-axis.
      const float t_at_depth = -depth / ray_direction_eye[2];
      const Point3f position_eye(Point3f::Zero() +
                                 t_at_depth * ray_direction_eye);
      return ion::math::ProjectPoint(world_from_eye_, position_eye);
    }
    case DepthType::kRayDepth: {
      DCHECK_GE(depth, 0.0f);
      DCHECK(std::isfinite(depth));
      const Vector3f normalized_ray_direction_eye = ion::math::Normalized(
          ion::math::ProjectPoint(eye_from_window_, Point3f(pixel, 0.0f)) -
          Point3f::Zero());
      const Point3f position_eye(Point3f::Zero() +
                                 normalized_ray_direction_eye * depth);
      return ion::math::ProjectPoint(world_from_eye_, position_eye);
    }
    default: {
      LOG(FATAL) << "Invalid or corrupt ProjectiveCamera. Ray end-point "
                    "computation is not supported with depth type "
                 << static_cast<int>(depth_type_) << ". "
                 << "Possible version mismatch.";
      return Point3f::Zero();
    }
  }
}

bool ProjectiveCamera::operator==(const ProjectiveCamera& other) const {
  return (other.image_size_ == image_size_ &&
          other.clip_from_eye_ == clip_from_eye_ &&
          other.eye_from_clip_ == eye_from_clip_ &&
          other.eye_from_world_ == eye_from_world_ &&
          other.world_from_eye_ == world_from_eye_ &&
          other.window_from_eye_ == window_from_eye_ &&
          other.eye_from_window_ == eye_from_window_ &&
          other.window_from_world_ == window_from_world_ &&
          other.world_from_window_ == world_from_window_ &&
          other.depth_type_ == depth_type_);
}

void ProjectiveCamera::ComputeNearFar(float* near_clip, float* far_clip) const {
  if (std::abs(clip_from_eye_[3][2]) < kOrthographicProjectionEpsilon) {
    // Orthographic projection
    *near_clip = (clip_from_eye_[2][3] + 1.0f) / clip_from_eye_[2][2];
    *far_clip = *near_clip - (2.0f / clip_from_eye_[2][2]);
  } else {
    // Perspective projection
    const Matrix4f& m = clip_from_eye_;

    // Test if the matrix is consistent with a finite projection.
    float finite_near_clip = (2.0f * m[2][3]) / (2.0f * m[2][2] - 2.0f);
    float finite_far_clip =
        ((m[2][2] - 1.0f) * finite_near_clip) / (m[2][2] + 1.0f);
    if (finite_near_clip > 0.0f && finite_far_clip > finite_near_clip &&
        finite_far_clip < std::numeric_limits<float>::max()) {
      *near_clip = finite_near_clip;
      *far_clip = finite_far_clip;
    } else {
      // If that fails, assume the matrix represents an infinite far clip
      // projection. Calculate the epsilon, if any, on the far clip plane.
      const float far_z_epsilon = m[2][2] + 1.0f;
      *near_clip = m[2][3] / (-2.0f + far_z_epsilon);
      *far_clip = std::numeric_limits<float>::max();
    }
  }
}

void ProjectiveCamera::ComputeClipPlanes(float* left, float* right,
                                         float* bottom, float* top,
                                         float* near_clip,
                                         float* far_clip) const {
  ComputeNearFar(near_clip, far_clip);
  const Matrix4f& m = clip_from_eye_;
  if (std::abs(clip_from_eye_[3][2]) < kOrthographicProjectionEpsilon) {
    // Orthographic projection
    *left = (-m[0][3] - 1.0f) / m[0][0];
    *right = (-m[0][3] + 1.0f) / m[0][0];
    *bottom = (-m[1][3] - 1.0f) / m[1][1];
    *top = (-m[1][3] + 1.0f) / m[1][1];
  } else {
    // Perspective projection
    *left = *near_clip * (m[0][2] - 1.0f) / m[0][0];
    *right = *near_clip * (m[0][2] + 1.0f) / m[0][0];
    *bottom = *near_clip * (m[1][2] - 1.0f) / m[1][1];
    *top = *near_clip * (m[1][2] + 1.0f) / m[1][1];
  }
}

}  // namespace base
}  // namespace seurat
