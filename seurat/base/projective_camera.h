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

#ifndef VR_SEURAT_BASE_PROJECTIVE_CAMERA_H_
#define VR_SEURAT_BASE_PROJECTIVE_CAMERA_H_

#include "ion/math/matrix.h"
#include "ion/math/matrixutils.h"
#include "ion/math/vector.h"
#include "seurat/base/camera.h"

namespace seurat {
namespace base {

// This class represents the intrinsic and extrinsic parameters of a projective
// implementation of a Camera.
class ProjectiveCamera : public Camera {
 public:
  // Defines various interpretations of depth values, used to compute ray end
  // points.
  //
  // Keep in sync with seurat::api::proto::DepthType (defined in
  // api/camera.proto).
  enum class DepthType {
    // Unspecified depth type.
    kUnspecified = 0,
    // Depths are the window-space Z coordinate (Z/W, as in Z buffer from GL) in
    // the range [0, 1].
    kWindowZ,
    // Depths are the negated eye-space Z coordinate in the range [0, inf).
    kEyeZ,
    // Depths are distances along a normalized ray (unit length direction
    // vector) through a pixel center. In other words, this is the distance
    // between the point and the origin in eye-space.
    kRayDepth
  };

  ProjectiveCamera() = default;

  // Creates a camera with the given image ("sensor") size in pixels, a
  // projection or clip space transformation in |clip_from_eye|, the view matrix
  // in |eye_from_world| and the corresponding inverse, or camera matrix, in
  // |world_from_eye|, and depth_type.
  ProjectiveCamera(const ion::math::Vector2i& image_size,
                   const ion::math::Matrix4f& clip_from_eye,
                   const ion::math::Matrix4f& eye_from_world,
                   const ion::math::Matrix4f& world_from_eye,
                   DepthType depth_type);

  // Creates a camera with the given image ("sensor") size in pixels,
  // |clip_from_eye|, |eye_from_world| and depth_type.
  ProjectiveCamera(const ion::math::Vector2i& image_size,
                   const ion::math::Matrix4f& clip_from_eye,
                   const ion::math::Matrix4f& eye_from_world,
                   DepthType depth_type)
      : ProjectiveCamera(image_size, clip_from_eye, eye_from_world,
                         ion::math::Inverse(eye_from_world), depth_type) {}

  // Creates a camera with the given image ("sensor") size in pixels,
  // |clip_from_eye|, |eye_from_world| and window-space Z depth interpretation.
  ProjectiveCamera(const ion::math::Vector2i& image_size,
                   const ion::math::Matrix4f& clip_from_eye,
                   const ion::math::Matrix4f& eye_from_world)
      : ProjectiveCamera(image_size, clip_from_eye, eye_from_world,
                         DepthType::kWindowZ) {}

  ~ProjectiveCamera() override = default;

  // Camera implementation.
  ion::math::Vector2i GetImageSize() const override;
  ion::math::Matrix4f GetWorldFromEye() const override;
  ion::math::Point3f RayOriginFloat(
      const ion::math::Point2f& pixel) const override;
  ion::math::Vector3f RayDirectionFloat(
      const ion::math::Point2f& pixel) const override;
  ion::math::Point3f RayEndFloat(const ion::math::Point2f& pixel,
                                 float depth) const override;

  // Returns true if the |other| camera is equal to this one.
  bool operator==(const ProjectiveCamera& other) const;

  // Returns true if the |other| camera is not equal to this one.
  bool operator!=(const ProjectiveCamera& other) const {
    return !(*this == other);
  }

  // Returns the clip-from-eye matrix, that transforms from eye coordinates to
  // clip coordinates.
  const ion::math::Matrix4f& GetClipFromEye() const { return clip_from_eye_; }

  // Returns the eye-from-clip matrix, that transforms from clip coordinates to
  // eye coordinates.
  const ion::math::Matrix4f& GetEyeFromClip() const { return eye_from_clip_; }

  // Returns the eye-from-world matrix, that transforms from world coordinates
  // to eye coordinates.
  const ion::math::Matrix4f& GetEyeFromWorld() const { return eye_from_world_; }

  // Returns the window-from-eye matrix, that transforms from eye coordinates
  // to window coordinates.
  const ion::math::Matrix4f& GetWindowFromEye() const {
    return window_from_eye_;
  }

  // Returns the eye-from-window matrix, that transforms from window coordinates
  // to eye coordinates.
  const ion::math::Matrix4f& GetEyeFromWindow() const {
    return eye_from_window_;
  }

  // Returns the window-from-world matrix, that transforms from world
  // coordinates to window coordinates.
  const ion::math::Matrix4f& GetWindowFromWorld() const {
    return window_from_world_;
  }

  // Returns the world-from-window matrix, that transforms from window
  // coordinates to world coordinates.
  const ion::math::Matrix4f& GetWorldFromWindow() const {
    return world_from_window_;
  }

  // Returns the depth type of the camera.
  DepthType GetDepthType() const { return depth_type_; }

  // Computes the near and far clipping planes from the clip-from-eye matrix.
  void ComputeNearFar(float* near_clip, float* far_clip) const;

  // Computes all six clip planes from the clip-from-eye matrix.
  void ComputeClipPlanes(float* left, float* right, float* top, float* bottom,
                         float* near_clip, float* far_clip) const;

 private:
  // Size of the image ("sensor") in pixels.
  ion::math::Vector2i image_size_;

  // The matrix that transforms eye coordinates to clip coordinates.
  ion::math::Matrix4f clip_from_eye_;

  // The matrix that transforms clip coordinates to eye coordinates.
  ion::math::Matrix4f eye_from_clip_;

  // The matrix that transforms world coordinates to eye coordinates.
  ion::math::Matrix4f eye_from_world_;

  // World-from-eye matrix (extrinsic camera parameters).
  ion::math::Matrix4f world_from_eye_;

  // The matrix that transforms eye coordinates to window coordinates.
  ion::math::Matrix4f window_from_eye_;

  // The matrix that transforms window coordinates to eye coordinates.
  ion::math::Matrix4f eye_from_window_;

  // The matrix that transforms from world coordinates to window coordinates.
  ion::math::Matrix4f window_from_world_;

  // The matrix that transforms from window coordinates to world coordinates.
  ion::math::Matrix4f world_from_window_;

  // The way this camera interprets depth values.
  DepthType depth_type_;
};

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_PROJECTIVE_CAMERA_H_
