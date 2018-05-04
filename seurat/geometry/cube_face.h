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

#ifndef VR_SEURAT_GEOMETRY_CUBE_FACE_H_
#define VR_SEURAT_GEOMETRY_CUBE_FACE_H_

#include <string>

#include "ion/math/matrix.h"
#include "ion/math/vector.h"

namespace seurat {
namespace geometry {

// The cube map coordinate space is equivalent to OpenGL camera space. Therefore
// the kFront face is directed toward the negative Z-axis.
enum class CubeFace : int { kFront = 0, kBack, kLeft, kRight, kBottom, kTop };

// Converts a string to its corresponding CubeFace enum value as follows:
//   "front" -> kFront
//   "back" -> kBack
//   "left" -> kLeft
//   "right" -> kRight
//   "bottom" -> kBottom
//   "top" -> kTop
CubeFace CubeFaceFromString(const std::string& name);

// Convert a cube-space location |position_cube_space| into the corresponding
// cube face value. Use, for example, to bucket point cloud samples to cube map
// faces.
CubeFace CubeFaceFromPosition(const ion::math::Point3f& position_cube_space);

// Returns whether the given string is a valid CubeFace string.
bool IsCubeFaceString(const std::string& name);

// Converts a CubeFace to a string.
//
// See CubeFaceFromString for more details.
const std::string& StringFromCubeFace(CubeFace face);

// Returns a 4x4 viewing matrix based on the given cube-face.
//
// This follows the convention of the old gluLookAt() function.  That is, the
// returned matrix transforms the given |face| in world-space (assuming "front"
// is facing the negative-z axis) to be in front of the camera in eye-space.
ion::math::Matrix4f LookAtMatrixFromFace(CubeFace face);

// Definition of cube-map coordinate system.
//
// UV coordinates are defined to be consistent with LookAtMatrixFromFace. The
// origin of the UV space is in the lower-left corner when looking through the
// cube-face camera. Positive U points to the right and positive V points up.
//
// Right  (positive X): U = +Z, V = +Y, Origin = ( 1, -1, -1)
// Left   (negative X): U = -Z, V = +Y, Origin = (-1, -1,  1)
// Top    (positive Y): U = +X, V = +Z, Origin = (-1,  1, -1)
// Bottom (negative Y): U = +X, V = -Z, Origin = (-1, -1,  1)
// Back   (positive Z): U = -X, V = +Y, Origin = ( 1, -1,  1)
// Front  (negative Z): U = +X, V = +Y, Origin = (-1, -1, -1)

// Returns the world-space position on the cube [-1,1]^3 for the given |uv|
// coordinate on the given |cube_face|.
ion::math::Point3f XYZFromUV(const ion::math::Point2f& uv, CubeFace cube_face);

// Returns the uv-coordinate on the given |cube_face| for the given world-space
// position |xyz|. |xyz| is assumed to be in the view frustum of the
// |cube_face|.
ion::math::Point2f UVFromXYZ(const ion::math::Point3f& xyz, CubeFace cube_face);

// Returns the world-space position for the given projective (homogeneous) |uvw|
// coordinate on the given |cube_face|. Note that this is a true world-space
// point, NOT a point on the [-1,1]^3 cube.
ion::math::Point3f XYZFromUVW(const ion::math::Point3f& uvw,
                              CubeFace cube_face);

// Returns the projective (homogenous) uvw coordinate on the given |cube_face|
// for the given world-space position |xyz|. |xyz| is assumed to be in the view
// frustum of the |cube_face|.
ion::math::Point3f UVWFromXYZ(const ion::math::Point3f& xyz,
                              CubeFace cube_face);

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_CUBE_FACE_H_
