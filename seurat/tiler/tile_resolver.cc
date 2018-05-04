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

#include "seurat/tiler/tile_resolver.h"

#include <numeric>

#include "seurat/base/ion_util_no_gl.h"
#include "seurat/geometry/convex_hull2d.h"
#include "seurat/geometry/plane.h"
#include "seurat/geometry/polygon.h"
#include "seurat/tiler/tile.h"

namespace seurat {
namespace tiler {

using geometry::Plane3f;
using geometry::Quad2f;
using geometry::Quad3f;
using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Range2f;
using ion::math::Vector2f;
using ion::math::Vector3f;

namespace {

// Returns true if all coordinates are finite values, false if any are Inf or
// NaN.
bool IsWellFormedQuad(const Quad3f& quad) {
  for (int i = 0; i < 4; ++i) {
    for (int d = 0; d < 3; ++d) {
      if (!std::isfinite(quad[i][d])) {
        return false;
      }
    }
  }
  return true;
}

}  // namespace

void RailTileResolver::Init(const PointSet& point_set) {
  subdivision_->Init(point_set);
}

bool RailTileResolver::Resolve(const BuildPartition& partition,
                               Tile* tile) const {
  std::array<Vector3f, 4> corner_rays =
      subdivision_->GetRails(partition.GetModel().cell);
  CHECK_GE(partition.GetModel().cell, 0);

  // Intersect these 3D rays with the partition plane to get the tile.
  const Plane3f plane(partition.GetModel().GetPlane());
  for (int i = 0; i < 4; ++i) {
    float t_hit;
    Vector3f ray_dir = corner_rays[i];
    CHECK_GT(ion::math::LengthSquared(ray_dir), 0.0f);
    if (!plane.IntersectRay(Point3f::Zero(), ray_dir, &t_hit)) {
      // This is an ill-formed tile.
      return false;
    }
    tile->quad[i] = ray_dir * t_hit + Point3f::Zero();
  }
  tile->cell = partition.GetModel().cell;
  return IsWellFormedQuad(tile->quad);
}

}  // namespace tiler
}  // namespace seurat
