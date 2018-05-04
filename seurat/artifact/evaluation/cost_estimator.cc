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

#include "seurat/artifact/evaluation/cost_estimator.h"

#include <cmath>
#include <random>

#include "ion/math/vector.h"
#include "seurat/geometry/fibonacci_sphere.h"

namespace seurat {
namespace artifact {

using geometry::Mesh;
using geometry::Raytracer;
using ion::math::Point3f;
using ion::math::Vector3f;

// static
CostEstimator CostEstimator::Build(const std::vector<SceneGeometry>& geometry) {
  // Opaque geometry is batched into a single Raytracer instance, with a map
  // indicating cost for each triangle.
  //
  // Each piece of transparent geometry is maintained separately.

  Mesh total_opaque_geometry(1);
  std::map<int, float> opaque_triangle_cost;

  std::vector<TransparentSceneGeometry> transparent;

  for (const auto& scene_geometry : geometry) {
    if (scene_geometry.opaque) {
      // Add to our batched opaque geometry.
      int base_index = total_opaque_geometry.GetTriangleCount();
      total_opaque_geometry.AppendMesh(scene_geometry.mesh);
      opaque_triangle_cost[base_index] = scene_geometry.cost;
    } else {
      TransparentSceneGeometry processed = {
          scene_geometry.cost, geometry::Raytracer::Build(scene_geometry.mesh)};
      transparent.push_back(std::move(processed));
    }
  }

  return CostEstimator(Raytracer::Build(total_opaque_geometry),
                       std::move(opaque_triangle_cost), std::move(transparent));
}

float CostEstimator::EstimateRayCost(const Point3f& ray_origin,
                                     const Vector3f& ray_direction) {
  // Estimate cost by first tracing a single ray through the batched opaque
  // geometry.
  //
  // If there is an intersection, then use the t_hit to limit the search through
  // each batch of transparent geometry.

  const int kMaxIntersectionCount = std::numeric_limits<int>::max();

  float t_hit_opaque = std::numeric_limits<float>::max();
  int triangle_index_hit_opaque;

  bool hit_opaque = opaque_->FindFirstHit(
      ray_origin, ray_direction, &t_hit_opaque, &triangle_index_hit_opaque);

  float cost_opaque = 0.0f;
  if (hit_opaque) {
    auto start_of_range =
        --opaque_triangle_cost_.upper_bound(triangle_index_hit_opaque);
    cost_opaque = start_of_range->second;
  }

  float cost_transparent = 0.0f;
  for (const auto& geometry : transparent_) {
    int intersections = geometry.raytracer->CountIntersections(
        ray_origin, ray_direction, t_hit_opaque, kMaxIntersectionCount);
    cost_transparent += geometry.cost * intersections;
  }

  return cost_opaque + cost_transparent;
}

float CostEstimator::EstimateCost(const ion::math::Range3f& headbox,
                                  int resolution, int samples_per_cell) {
  std::mt19937 random;
  std::uniform_real_distribution<float> unit(0.0f, 1.0f);

  const Point3f headbox_origin = headbox.GetMinPoint();
  const Vector3f headbox_size = headbox.GetSize();

  float total_cost = 0.0f;

  for (int x = 0; x < resolution; ++x) {
    for (int y = 0; y < resolution; ++y) {
      for (int z = 0; z < resolution; ++z) {
        float fibonacci_sphere_scrambler = unit(random) * 2.0f * M_PI;
        for (int i = 0; i < samples_per_cell; ++i) {
          // Offset within the headbox in the range [0, resolution]^3
          Vector3f offset = Vector3f(x, y, z) +
                            Vector3f(unit(random), unit(random), unit(random));

          Point3f ray_origin =
              headbox_origin + offset * headbox_size / resolution;

          // Sample using fibonacci spiral gridpoints on a sphere.
          Vector3f ray_direction =
              Point3f(geometry::GenerateFibonacciSpherePoint(
                  samples_per_cell, fibonacci_sphere_scrambler, i)) -
              Point3f::Zero();

          total_cost += EstimateRayCost(ray_origin, ray_direction);
        }
      }
    }
  }

  int total_samples = resolution * resolution * resolution * samples_per_cell;
  return total_cost / total_samples;
}

}  // namespace artifact
}  // namespace seurat
