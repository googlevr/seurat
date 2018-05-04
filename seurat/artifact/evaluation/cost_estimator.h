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

#ifndef VR_SEURAT_ARTIFACT_EVALUATION_COST_ESTIMATOR_H_
#define VR_SEURAT_ARTIFACT_EVALUATION_COST_ESTIMATOR_H_

#include <map>

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "seurat/geometry/mesh.h"
#include "seurat/geometry/raytracer.h"

namespace seurat {
namespace artifact {

// Estimates the cost of rendering a scene consisting of opaque & non-opaque
// triangle-meshes, each with their own cost-factor.
//
// The cost model returns the estimated integral over all rays from a headbox of
// the cost of that ray, computed as the sum of the cost of geometry intersected
// until an opaque surface is hit.
class CostEstimator {
 public:
  struct SceneGeometry {
    // A triangle mesh.
    geometry::Mesh mesh;

    // A measure of the cost of shaded fragments of this scene geometry.
    float cost;

    // Whether this scene geometry is opaque or transparent.
    bool opaque;
  };

  // Builds a CostEstimator for the given scene geometry.
  static CostEstimator Build(const std::vector<SceneGeometry>& geometry);

  // The cost of all fragments which would be shaded to process the ray from
  // |ray_origin| in |ray_direction|.
  float EstimateRayCost(const ion::math::Point3f& ray_origin,
                        const ion::math::Vector3f& ray_direction);

  // The mean cost for rays randomly generated within the |headbox|.
  //
  // |samples_per_cell| are cast randomly from ray origins generated within
  // uniformly-spaced cells within the |headbox| at the specified |resolution|.
  //
  // The total number of rays cast is resolution^3 * samples_per_cell.
  float EstimateCost(const ion::math::Range3f& headbox, int resolution,
                     int samples_per_cell);

 private:
  // Internal representation of scene geometry, for transparent geometry.
  struct TransparentSceneGeometry {
    // Cost of rasterizing any fragment of the scene of the 'raytracer'.
    float cost;
    // Traces rays through scene geometry.
    std::unique_ptr<geometry::Raytracer> raytracer;
  };

  CostEstimator(std::unique_ptr<geometry::Raytracer> opaque,
                std::map<int, float> opaque_triangle_cost,
                std::vector<TransparentSceneGeometry> transparent)
      : opaque_(std::move(opaque)),
        opaque_triangle_cost_(std::move(opaque_triangle_cost)),
        transparent_(std::move(transparent)) {}

  // A single scene for all opaque geometry.
  std::unique_ptr<geometry::Raytracer> opaque_;

  // 'opaque_' traces rays into a scene of triangles, each with an index.
  //
  // Consecutive ranges of these triangles have identical cost.
  //
  // This stores a mapping from the first index of each range to the cost of
  // each triangle in that range.
  std::map<int, float> opaque_triangle_cost_;

  // A raytracing scene for each set of transparent geometry.
  std::vector<TransparentSceneGeometry> transparent_;
};

}  // namespace artifact
}  // namespace seurat

#endif  // VR_SEURAT_ARTIFACT_EVALUATION_COST_ESTIMATOR_H_
