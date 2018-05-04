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

#include "seurat/ingest/clamping_view_group_loader.h"

#include "ion/math/vector.h"
#include "seurat/base/camera.h"
#include "seurat/base/status.h"
#include "seurat/image/ldi.h"

namespace seurat {
namespace ingest {

using base::Camera;
using image::Ldi4f;
using ion::math::Point3f;
using ion::math::Vector3f;

namespace {

// Wraps an arbitrary Camera with unified skybox clamping behavior.
//
// Points beyond the skybox are clamped to lie on the skybox.
class ClampingCamera : public base::Camera {
 public:
  ClampingCamera(float skybox_radius, bool zero_is_infinite,
                 std::shared_ptr<base::Camera> delegate)
      : skybox_radius_(skybox_radius),
        zero_is_infinite_(zero_is_infinite),
        delegate_(std::move(delegate)) {}
  ~ClampingCamera() override = default;

  // Camera implementation.
  ion::math::Vector2i GetImageSize() const override {
    return delegate_->GetImageSize();
  }
  ion::math::Matrix4f GetWorldFromEye() const override {
    return delegate_->GetWorldFromEye();
  }
  ion::math::Point3f RayOriginFloat(
      const ion::math::Point2f& pixel) const override {
    return delegate_->RayOriginFloat(pixel);
  }
  ion::math::Vector3f RayDirectionFloat(
      const ion::math::Point2f& pixel) const override {
    return delegate_->RayDirectionFloat(pixel);
  }

  // This is where the real action happens.
  //
  // Invokes the wrapped camera & clamps the ray end.
  ion::math::Point3f RayEndFloat(const ion::math::Point2f& pixel,
                                 float depth) const override;

 private:
  // Half the side length of an origin-centered cube defining the sky box.
  const float skybox_radius_;

  // If true, depths of zero will be considered infinite.
  //
  // This is useful, for example, when processing output from a game engine
  // which cannot properly represent points at infinity in its depth buffer.
  const bool zero_is_infinite_;

  // The camera to wrap.
  const std::shared_ptr<base::Camera> delegate_;
};

ion::math::Point3f ClampingCamera::RayEndFloat(const ion::math::Point2f& pixel,
                                               float depth) const {
  Point3f original_point = delegate_->RayEndFloat(pixel, depth);
  bool use_original = true;
  if (zero_is_infinite_ && depth == 0.0f) {
    use_original = false;
  }
  for (int d = 0; d < 3; ++d) {
    if (std::fabs(original_point[d]) > skybox_radius_) {
      use_original = false;
    }
  }
  if (use_original) {
    // This is the common case.
    return original_point;
  }

  Point3f origin = RayOriginFloat(pixel);
  Vector3f direction = RayDirectionFloat(pixel);

  float closest_positive_scaler = std::numeric_limits<float>::infinity();
  for (int d = 0; d < 3; ++d) {
    // Solve for the scaler which results in the ray intersecting one of the
    // skybox planes.
    //
    // In other words, we want to solve for hit_scalers:
    //   abs(origin[d] + direction[d] * hit_scalers[d]) == skybox_radius
    float scaler;
    scaler = (skybox_radius_ - origin[d]) / direction[d];
    if (std::isfinite(scaler) && scaler > 0.0f) {
      closest_positive_scaler = std::min(scaler, closest_positive_scaler);
    }
    scaler = (-skybox_radius_ - origin[d]) / direction[d];
    if (std::isfinite(scaler) && scaler > 0.0f) {
      closest_positive_scaler = std::min(scaler, closest_positive_scaler);
    }
  }

  return origin + direction * closest_positive_scaler;
}

}  // namespace

base::Status ClampingViewGroupLoader::LoadViewGroup(
    int view_group_index, std::vector<std::shared_ptr<Camera>>* cameras,
    std::vector<Ldi4f>* ldis) const {
  base::Status status =
      delegate_->LoadViewGroup(view_group_index, cameras, ldis);

  if (cameras != nullptr) {
    for (auto& camera_ptr : *cameras) {
      camera_ptr = WrapCamera(std::move(camera_ptr));
    }
  }

  return status;
}

std::shared_ptr<base::Camera> ClampingViewGroupLoader::WrapCamera(
    std::shared_ptr<base::Camera> original) const {
  return std::make_shared<ClampingCamera>(skybox_radius_, zero_is_infinite_,
                                          std::move(original));
}

}  // namespace ingest
}  // namespace seurat
