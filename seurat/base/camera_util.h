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

#ifndef VR_SEURAT_BASE_CAMERA_UTIL_H_
#define VR_SEURAT_BASE_CAMERA_UTIL_H_

#include <memory>

#include "seurat/api/camera.pb.h"
#include "seurat/base/projective_camera.h"

namespace seurat {
namespace base {

// Returns a ProjectiveCamera read from protocol buffer object |proto|.
std::unique_ptr<Camera> CameraFromProto(const api::proto::Camera& proto);

// Wraps another camera with the specified |translation|.
std::unique_ptr<Camera> TranslateCamera(const ion::math::Vector3f& translation,
                                        std::shared_ptr<base::Camera> original);

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_CAMERA_UTIL_H_
