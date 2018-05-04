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

#include "seurat/base/math_util.h"

#include "absl/strings/str_cat.h"

namespace seurat {
namespace base {

ion::math::Vector2i Vector2iFromProto(const api::proto::Vector2i& proto) {
  return ion::math::Vector2i(proto.x(), proto.y());
}

void Vector2iToProto(const ion::math::Vector2i& v,
                     api::proto::Vector2i* proto) {
  proto->set_x(v[0]);
  proto->set_y(v[1]);
}

ion::math::Point3f Point3fFromProto(const api::proto::Point3f& proto) {
  return ion::math::Point3f(proto.x(), proto.y(), proto.z());
}

void Point3fToProto(const ion::math::Point3f& v, api::proto::Point3f* proto) {
  proto->set_x(v[0]);
  proto->set_y(v[1]);
  proto->set_z(v[2]);
}

ion::math::Matrix4f Matrix4fFromProto(const api::proto::Matrix4f& proto) {
  CHECK_EQ(16, proto.m_size())
      << "Matrix4f protobuf must have exactly 16 elements.";

  return ion::math::Matrix4f(
      proto.m(0), proto.m(1), proto.m(2), proto.m(3),    //
      proto.m(4), proto.m(5), proto.m(6), proto.m(7),    //
      proto.m(8), proto.m(9), proto.m(10), proto.m(11),  //
      proto.m(12), proto.m(13), proto.m(14), proto.m(15));
}

void Matrix4fToProto(const ion::math::Matrix4f& m,
                     api::proto::Matrix4f* proto) {
  for (int row = 0; row < 4; ++row) {
    for (int col = 0; col < 4; ++col) {
      proto->add_m(m[row][col]);
    }
  }
}

base::Status ValidateMatrix4f(const Json::Value& matrix) {
  if (!matrix.isArray()) {
    return base::InvalidArgumentError(
        absl::StrCat("matrix is not an array: ", matrix.toStyledString()));
  }

  if (matrix.size() != 16) {
    return base::InvalidArgumentError(
        absl::StrCat("size of matrix != 16: ", matrix.toStyledString()));
  }

  for (int i = 0; i < 16; ++i) {
    if (!matrix[i].isNumeric()) {
      return base::InvalidArgumentError(absl::StrCat(
          "matrix[", i, "] is not numeric: ", matrix[i].toStyledString()));
    }
  }

  return base::OkStatus();
}

Json::Value JsonFromMatrix4f(const ion::math::Matrix4f& matrix) {
  Json::Value value(Json::arrayValue);
  for (int row = 0; row < 4; ++row) {
    for (int col = 0; col < 4; ++col) {
      value.append(matrix(row, col));
    }
  }
  return value;
}

ion::math::Matrix4f Matrix4fFromJson(const Json::Value& matrix) {
  DCHECK(ValidateMatrix4f(matrix).ok());
  float m[16];
  for (int i = 0; i < 16; ++i) {
    m[i] = matrix[i].asFloat();
  }
  return ion::math::Matrix4f(m[0], m[1], m[2], m[3],    //
                             m[4], m[5], m[6], m[7],    //
                             m[8], m[9], m[10], m[11],  //
                             m[12], m[13], m[14], m[15]);
}

}  // namespace base
}  // namespace seurat
