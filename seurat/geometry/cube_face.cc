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

#include "seurat/geometry/cube_face.h"

#include <algorithm>
#include <map>

#include "ion/math/matrixutils.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "seurat/base/ion_util_no_gl.h"

namespace seurat {
namespace geometry {

using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Point4f;
using ion::math::Vector3f;

namespace {

const std::map<std::string, CubeFace>& GetStringFaceMap() {
  static const auto m = new std::map<std::string, CubeFace>{
      {"front", CubeFace::kFront},   {"back", CubeFace::kBack},
      {"left", CubeFace::kLeft},     {"right", CubeFace::kRight},
      {"bottom", CubeFace::kBottom}, {"top", CubeFace::kTop}};
  return *m;
}

const std::map<CubeFace, std::string>& GetFaceStringMap() {
  static const auto m = [] {
    auto rev = new std::map<CubeFace, std::string>();
    for (const auto& kv : GetStringFaceMap()) {
      rev->insert({kv.second, kv.first});
    }
    return rev;
  }();
  return *m;
}

}  // namespace

CubeFace CubeFaceFromString(const std::string& name) {
  const auto& faces = GetStringFaceMap();
  auto iter = faces.find(name);
  CHECK(iter != faces.end()) << "Invalid CubeFace name: " << name;
  return iter->second;
}

CubeFace CubeFaceFromPosition(const Point3f& position_cube_space) {
  int ordered_signed_axis = base::MajorAxisFromPosition(position_cube_space);
  // Shift negative faces to give this total order: X, Y, Z, -X, -Y, -Z.
  if (position_cube_space[ordered_signed_axis] < 0.0f) {
    ordered_signed_axis += 3;
  }
  const CubeFace kFacesFromOrderedSignedAxes[6] = {
      CubeFace::kRight, CubeFace::kTop,    CubeFace::kBack,
      CubeFace::kLeft,  CubeFace::kBottom, CubeFace::kFront};
  return kFacesFromOrderedSignedAxes[ordered_signed_axis];
}

bool IsCubeFaceString(const std::string& name) {
  return GetStringFaceMap().count(name) > 0;
}

const std::string& StringFromCubeFace(CubeFace face) {
  const auto& names = GetFaceStringMap();
  auto iter = names.find(face);
  DCHECK(iter != names.end());
  return iter->second;
}

ion::math::Matrix4f LookAtMatrixFromFace(CubeFace face) {
  switch (face) {
    // clang-format off
    case CubeFace::kFront:
      return Matrix4f( 1.0f,  0.0f,  0.0f,  0.0f,
                       0.0f,  1.0f,  0.0f,  0.0f,
                       0.0f,  0.0f,  1.0f,  0.0f,
                       0.0f,  0.0f,  0.0f,  1.0f);
    case CubeFace::kBack:
      return Matrix4f(-1.0f,  0.0f,  0.0f,  0.0f,
                       0.0f,  1.0f,  0.0f,  0.0f,
                       0.0f,  0.0f, -1.0f,  0.0f,
                       0.0f,  0.0f,  0.0f,  1.0f);
    case CubeFace::kLeft:
      return Matrix4f( 0.0f,  0.0f, -1.0f,  0.0f,
                       0.0f,  1.0f,  0.0f,  0.0f,
                       1.0f,  0.0f,  0.0f,  0.0f,
                       0.0f,  0.0f,  0.0f,  1.0f);
    case CubeFace::kRight:
      return Matrix4f( 0.0f,  0.0f,  1.0f,  0.0f,
                       0.0f,  1.0f,  0.0f,  0.0f,
                      -1.0f,  0.0f,  0.0f,  0.0f,
                       0.0f,  0.0f,  0.0f,  1.0f);
    case CubeFace::kBottom:
      return Matrix4f( 1.0f,  0.0f,  0.0f,  0.0f,
                       0.0f,  0.0f, -1.0f,  0.0f,
                       0.0f,  1.0f,  0.0f,  0.0f,
                       0.0f,  0.0f,  0.0f,  1.0f);
    case CubeFace::kTop:
      return Matrix4f( 1.0f,  0.0f,  0.0f,  0.0f,
                       0.0f,  0.0f,  1.0f,  0.0f,
                       0.0f, -1.0f,  0.0f,  0.0f,
                       0.0f,  0.0f,  0.0f,  1.0f);
    // clang-format on
    default:
      LOG(FATAL) << "Invalid value in switch statement: "
                 << static_cast<int>(face);
      return Matrix4f::Identity();
  }
}

Point3f XYZFromUV(const Point2f& uv, CubeFace cube_face) {
  // Convert range from [0,1] to [-1,1].
  const float u = 2.0f * uv[0] - 1.0f;
  const float v = 2.0f * uv[1] - 1.0f;

  switch (cube_face) {
    case CubeFace::kRight:
      return Point3f(1.0f, v, u);
    case CubeFace::kLeft:
      return Point3f(-1.0f, v, -u);
    case CubeFace::kTop:
      return Point3f(u, 1.0f, v);
    case CubeFace::kBottom:
      return Point3f(u, -1.0f, -v);
    case CubeFace::kBack:
      return Point3f(-u, v, 1.0f);
    case CubeFace::kFront:
      return Point3f(u, v, -1.0f);
  }
}

Point2f UVFromXYZ(const Point3f& xyz, CubeFace cube_face) {
  switch (cube_face) {
    case CubeFace::kRight: {
      const float inv_w = 1.0f / xyz[0];
      return Point2f(0.5f * (xyz[2] * inv_w + 1.0f),
                     0.5f * (xyz[1] * inv_w + 1.0f));
    }
    case CubeFace::kLeft: {
      const float inv_w = 1.0f / -xyz[0];
      return Point2f(0.5f * (-xyz[2] * inv_w + 1.0f),
                     0.5f * (xyz[1] * inv_w + 1.0f));
    }
    case CubeFace::kTop: {
      const float inv_w = 1.0f / xyz[1];
      return Point2f(0.5f * (xyz[0] * inv_w + 1.0f),
                     0.5f * (xyz[2] * inv_w + 1.0f));
    }
    case CubeFace::kBottom: {
      const float inv_w = 1.0f / -xyz[1];
      return Point2f(0.5f * (xyz[0] * inv_w + 1.0f),
                     0.5f * (-xyz[2] * inv_w + 1.0f));
    }
    case CubeFace::kBack: {
      const float inv_w = 1.0f / xyz[2];
      return Point2f(0.5f * (-xyz[0] * inv_w + 1.0f),
                     0.5f * (xyz[1] * inv_w + 1.0f));
    }
    case CubeFace::kFront: {
      const float inv_w = 1.0f / -xyz[2];
      return Point2f(0.5f * (xyz[0] * inv_w + 1.0f),
                     0.5f * (xyz[1] * inv_w + 1.0f));
    }
  }
}

Point3f XYZFromUVW(const Point3f& uvw, CubeFace cube_face) {
  // Convert range from [0,1] to [-1,1].
  const float w = uvw[2];
  const float u = 2.0f * uvw[0] - w;
  const float v = 2.0f * uvw[1] - w;

  switch (cube_face) {
    case CubeFace::kRight:
      return Point3f(w, v, u);
    case CubeFace::kLeft:
      return Point3f(-w, v, -u);
    case CubeFace::kTop:
      return Point3f(u, w, v);
    case CubeFace::kBottom:
      return Point3f(u, -w, -v);
    case CubeFace::kBack:
      return Point3f(-u, v, w);
    case CubeFace::kFront:
      return Point3f(u, v, -w);
  }
}

Point3f UVWFromXYZ(const Point3f& xyz, CubeFace cube_face) {
  switch (cube_face) {
    case CubeFace::kRight: {
      const float w = xyz[0];
      return Point3f(0.5f * (xyz[2] + w), 0.5f * (xyz[1] + w), w);
    }
    case CubeFace::kLeft: {
      const float w = -xyz[0];
      return Point3f(0.5f * (-xyz[2] + w), 0.5f * (xyz[1] + w), w);
    }
    case CubeFace::kTop: {
      const float w = xyz[1];
      return Point3f(0.5f * (xyz[0] + w), 0.5f * (xyz[2] + w), w);
    }
    case CubeFace::kBottom: {
      const float w = -xyz[1];
      return Point3f(0.5f * (xyz[0] + w), 0.5f * (-xyz[2] + w), w);
    }
    case CubeFace::kBack: {
      const float w = xyz[2];
      return Point3f(0.5f * (-xyz[0] + w), 0.5f * (xyz[1] + w), w);
    }
    case CubeFace::kFront: {
      const float w = -xyz[2];
      return Point3f(0.5f * (xyz[0] + w), 0.5f * (xyz[1] + w), w);
    }
  }
}

}  // namespace geometry
}  // namespace seurat
