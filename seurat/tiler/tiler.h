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

#ifndef VR_SEURAT_TILER_TILER_H_
#define VR_SEURAT_TILER_TILER_H_

#include <memory>
#include <vector>

#include "ion/math/vector.h"
#include "seurat/geometry/plane.h"
#include "seurat/tiler/point_set.h"
#include "seurat/tiler/tile.h"

namespace seurat {
namespace tiler {

// Generates tiles to approximate the geometry sampled by a PointSet.
class Tiler {
 public:
  virtual ~Tiler() = default;

  virtual std::vector<Tile> Run(const PointSet& point_set) = 0;
};

// Builds a Tiler.
class TilerFactory {
 public:
  struct Parameters {
    // The maximum number of tiles to generate.
    int tile_count;

    // The maximum amount of average overdraw.
    float overdraw_factor;

    // The maximum amount of peak overdraw, bounded by sampling over many
    // directions.
    float peak_overdraw_factor = 1000.0f;

    // The field-of-view of the cameras to use to sample peak overdraw.
    float peak_overdraw_field_of_view_degrees = 90.0f;

    // The number of camera poses to sample when bounding peak overdraw.
    //
    // If zero, peak overdraw is unbounded.
    int peak_overdraw_samples = 100;

    // The radius of the viewing volume.
    float headbox_radius = 0.5f;

    // The half-side-length of the skybox.
    float skybox_radius = 200.0f;

    // The maximum number of threads to use.
    int thread_count;

    // If true, prioritize speed over quality.
    bool fast = false;

    // The minimum amount to tessellate (i.e. the min depth of a quadtree).
    //
    // The resulting tile count will be at least
    //   4^[min_subdivision_level] * 6
    // for a complete scene.
    int min_subdivision_level = 3;

    // The maximum amount to tessellate (i.e. the maximum depth of a quadtree).
    int max_subdivision_level = 7;
  };

  // Returns an instance of the current-best tiler.
  static std::unique_ptr<Tiler> CreateDefaultTiler(
      const Parameters& parameters) {
    return CreateSelectionTiler(parameters);
  }

  // Returns an instance of a SelectionTiler.
  static std::unique_ptr<Tiler> CreateSelectionTiler(
      const Parameters& parameters);
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_TILER_H_
