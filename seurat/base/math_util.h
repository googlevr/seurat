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

#ifndef VR_SEURAT_BASE_MATH_UTIL_H_
#define VR_SEURAT_BASE_MATH_UTIL_H_

#include "ion/math/matrix.h"
#include "ion/math/vector.h"
#include "json/json.h"
#include "seurat/api/math.pb.h"
#include "seurat/base/status.h"

namespace seurat {
namespace base {

// Returns a Vector2i read from protocol buffer object |proto|.
ion::math::Vector2i Vector2iFromProto(const api::proto::Vector2i& proto);

// Writes |v| to protocol buffer object |proto|.
void Vector2iToProto(const ion::math::Vector2i& v, api::proto::Vector2i* proto);

// Returns a Point3f read from protocol buffer object |proto|.
ion::math::Point3f Point3fFromProto(const api::proto::Point3f& proto);

// Writes |p| to protocol buffer object |proto|.
void Point3fToProto(const ion::math::Point3f& p, api::proto::Point3f* proto);

// Returns a Matrix4f read from protocol buffer object |proto|.
ion::math::Matrix4f Matrix4fFromProto(const api::proto::Matrix4f& proto);

// Writes |m| to protocol buffer object |proto|.
void Matrix4fToProto(const ion::math::Matrix4f& m, api::proto::Matrix4f* proto);

// Validates if |matrix| is a correct JSON serialization of a Matrix4f.
base::Status ValidateMatrix4f(const Json::Value& matrix);

// Returns a JSON serialization of the |matrix|.
Json::Value JsonFromMatrix4f(const ion::math::Matrix4f& matrix);

// Deserializes |matrix| into a Matrix4f and returns it. Assumes that |matrix|
// is a valid serialization of a Matrix4f.
ion::math::Matrix4f Matrix4fFromJson(const Json::Value& matrix);

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_MATH_UTIL_H_
