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

#ifndef VR_SEURAT_TILER_GEOMETRY_MODEL_H_
#define VR_SEURAT_TILER_GEOMETRY_MODEL_H_

#include "ion/math/vector.h"
#include "seurat/geometry/plane.h"

namespace seurat {
namespace tiler {

// Parameters for the proxy geometry used to represent a set of points.
struct GeometryModel {
  GeometryModel()
      : cell(-1),
        center(ion::math::Point3f::Zero()),
        normal(ion::math::Vector3f::AxisZ()) {}

  // Returns a planar representation of this surface proxy.
  geometry::Plane3f GetPlane() const { return {center, normal}; }

  bool operator!=(const GeometryModel& rhs) const { return !(*this == rhs); }
  bool operator==(const GeometryModel& rhs) const {
    return cell == rhs.cell && center == rhs.center && normal == rhs.normal;
  }

  // The Subdivision cell in which this GeometryModel lives.
  int cell;

  // The "center" of this piece of geometry.  The precise definition of which is
  // left to the GeometrySolver.
  ion::math::Point3f center;

  // The normal to the plane for this piece of proxy geometry.
  ion::math::Vector3f normal;
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_GEOMETRY_MODEL_H_
