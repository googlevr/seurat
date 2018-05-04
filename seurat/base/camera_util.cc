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

#include "seurat/base/camera_util.h"

#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "seurat/api/camera.pb.h"
#include "seurat/base/projective_camera_util.h"

namespace seurat {
namespace base {

using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Vector2i;
using ion::math::Vector3f;

std::unique_ptr<Camera> CameraFromProto(const api::proto::Camera& proto) {
  if (proto.has_projective()) {
    return std::unique_ptr<Camera>(
        new ProjectiveCamera(ProjectiveCameraFromProto(proto.projective())));
  } else {
    LOG(FATAL) << "Unknown camera type in proto.";
    return nullptr;
  }
}

namespace {

// Wraps an arbitrary Camera with a translation.
class TranslatingCamera : public base::Camera {
 public:
  TranslatingCamera(Vector3f translation,
                    std::shared_ptr<base::Camera> delegate)
      : translation_(translation), delegate_(std::move(delegate)) {}
  ~TranslatingCamera() override = default;

  // Camera implementation.
  Vector2i GetImageSize() const override { return delegate_->GetImageSize(); }
  Matrix4f GetWorldFromEye() const override {
    return ion::math::TranslationMatrix(translation_) *
           delegate_->GetWorldFromEye();
  }
  Point3f RayOriginFloat(const Point2f& pixel) const override {
    return delegate_->RayOriginFloat(pixel) + translation_;
  }
  Vector3f RayDirectionFloat(const Point2f& pixel) const override {
    return delegate_->RayDirectionFloat(pixel);
  }
  Point3f RayEndFloat(const Point2f& pixel, float depth) const override {
    return delegate_->RayEndFloat(pixel, depth) + translation_;
  }

 private:
  // Translation vector.
  const Vector3f translation_;

  // The camera to wrap.
  const std::shared_ptr<base::Camera> delegate_;
};

}  // namespace

std::unique_ptr<Camera> TranslateCamera(
    const Vector3f& translation, std::shared_ptr<base::Camera> original) {
  return std::unique_ptr<Camera>(
      new TranslatingCamera(translation, std::move(original)));
}

}  // namespace base
}  // namespace seurat
