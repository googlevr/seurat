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

#include "ion/math/matrix.h"
#include "ion/math/vector.h"
#include "gtest/gtest.h"

namespace seurat {
namespace base {

using ion::math::Matrix4f;
using ion::math::Point3f;
using ion::math::Vector2i;

namespace {

TEST(MathUtil, Vector2iSerializeAndDeserialize) {
  const Vector2i v(3, 5);
  api::proto::Vector2i proto;
  Vector2iToProto(v, &proto);
  std::string proto_string = proto.SerializeAsString();
  api::proto::Vector2i deserialized_proto;
  deserialized_proto.ParseFromString(proto_string);
  Vector2i deserialized = Vector2iFromProto(deserialized_proto);
  EXPECT_EQ(v, deserialized);
}

TEST(MathUtil, Point3fSerializeAndDeserialize) {
  const Point3f v(3.0f, 5.0f, 7.0f);
  api::proto::Point3f proto;
  Point3fToProto(v, &proto);
  std::string proto_string = proto.SerializeAsString();
  api::proto::Point3f deserialized_proto;
  deserialized_proto.ParseFromString(proto_string);
  Point3f deserialized = Point3fFromProto(deserialized_proto);
  EXPECT_EQ(v, deserialized);
}

TEST(MathUtil, Matrix4fSerializeAndDeserializeProto) {
  const Matrix4f source_matrix(1.0f, 2.0f, 3.0f, 4.0f,     //
                               5.0f, 6.0f, 7.0f, 8.0f,     //
                               9.0f, 10.0f, 11.0f, 12.0f,  //
                               13.0f, 14.0f, 15.0f, 16.0f);
  api::proto::Matrix4f source_matrix_proto;
  Matrix4fToProto(source_matrix, &source_matrix_proto);
  std::string proto_string = source_matrix_proto.SerializeAsString();
  api::proto::Matrix4f deserialized_matrix_proto;
  deserialized_matrix_proto.ParseFromString(proto_string);
  Matrix4f deserialized_matrix = Matrix4fFromProto(deserialized_matrix_proto);
  EXPECT_EQ(source_matrix, deserialized_matrix);
}

TEST(MathUtil, Matrix4fSerializeAndDeserializeJson) {
  const Matrix4f matrix(1.0f, 2.0f, 3.0f, 4.0f,     //
                        5.0f, 6.0f, 7.0f, 8.0f,     //
                        9.0f, 10.0f, 11.0f, 12.0f,  //
                        13.0f, 14.0f, 15.0f, 16.0f);
  Json::Value matrix_json = JsonFromMatrix4f(matrix);
  EXPECT_TRUE(ValidateMatrix4f(matrix_json).ok());
  const Matrix4f matrix_copy = Matrix4fFromJson(matrix_json);
  EXPECT_EQ(matrix, matrix_copy);
}

}  // namespace
}  // namespace base
}  // namespace seurat
