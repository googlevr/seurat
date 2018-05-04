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

#ifndef VR_SEURAT_BAKER_FRAMEWORK_FRAME_H_
#define VR_SEURAT_BAKER_FRAMEWORK_FRAME_H_

#include "ion/math/matrix.h"
#include "absl/types/span.h"
#include "seurat/base/array2d.h"
#include "seurat/base/color.h"
#include "seurat/geometry/plane.h"
#include "seurat/geometry/quad.h"

namespace seurat {
namespace baker {

// A quad to be rendered at runtime in a particular order.
//
// A collection of Frames forms a polygon soup to be used as coarse
// proxy-geometry for image-based rendering.
//
// For example, a collection of Frames may have textures assigned to them based
// on projective texture mapping, where samples are projected towards the origin
// of the world, which is assumed to be the primary location from which the
// scene is viewed.
struct Frame {
  // The order in which this quad will be rendered (0 indicates the first frame
  // to be rendered).
  int draw_order;

  // The geometry of the frame.
  geometry::Quad3f quad;

  // The w-component of the texture coordinates for each quad vertex.
  std::array<float, 4> texcoord_w;
};

bool operator==(const Frame& a, const Frame& b);
inline bool operator!=(const Frame& a, const Frame& b) { return !(a == b); }

// Initializes (or resets) the draw_order of all frames according to
// approximate back-to-front sorting, relative to the origin.
void InitializeApproximateDrawOrder(absl::Span<Frame> frames);

// Returns the plane of the quad for a given |frame|.
geometry::Plane3f PlaneFromFrame(const Frame& frame);

// Transforms the given frame into a coplanar frame which is dilated by an
// amount corresponding to approximately 1 pixel when viewed from the origin.
bool DilateFrame(float resolution, Frame* frame);

// World-space to frame-space.
//
// Transforms a 3D point on a frame onto normalized 2D coordinates.
//
// Frame coordinates are homogeneous (projective) texture coordinates that are
// in the range [0, 1]x[0, 1] after dividing by W.
//
// (0, w3, w3)  (w2, w2, w2)
//        3-------2
//        |     / |
//        |    /  |
//        |   /   |
//        |  /    |
//        | /     |
//        0-------1
// (0, 0, w0)  (w1, 0, w1)
//
// Note: Results are *undefined* if the point lies outside the frame.
// Furthermore, due to floating-point rounding, it is possible that edge-points
// result in normalized frame coordinates *outside* the unit-square.
//
// Returns whether |point_world| is inside the frame.
bool WorldToFrame(const Frame& frame, const ion::math::Point3f& point_world,
                  ion::math::Point3f* point_frame);

// Frame-space to world-space.
//
// This is the inverse of WorldToFrame.
//
// Returns whether |point_frame| is inside the frame.
bool FrameToWorld(const Frame& frame, const ion::math::Point3f& point_frame,
                  ion::math::Point3f* point_world);

// Converts a "freespace" ray defined by the given |start| and |direction|
// into its corresponding frame-space point.
//
// Returns false if there is no corresponding frame_space point.
//
// This is the intersection of the ray from |start| towards |direction| with the
// |frame|.
bool FreespaceRayToFrameSpace(const Frame& frame,
                              const ion::math::Point3f& start,
                              const ion::math::Vector3f& direction,
                              ion::math::Point2f* frame_space);

// Converts a "solid" ray defined by the given |start| and |end| points
// into its corresponding frame-space point.
//
// Returns false if there is no corresponding frame_space point.
//
// This is the intersection of the ray from (0,0,0)->|end| with the |frame|.
//
// Note that this is consistent with the RayClassifier.
bool SolidRayToFrameSpace(const Frame& frame, const ion::math::Point3f& start,
                          const ion::math::Point3f& end,
                          ion::math::Point2f* frame_space);

}  // namespace baker
}  // namespace seurat

#endif  // VR_SEURAT_BAKER_FRAMEWORK_FRAME_H_
