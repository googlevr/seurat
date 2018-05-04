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

#ifndef VR_SEURAT_GEOMETRY_RAYTRACER_H_
#define VR_SEURAT_GEOMETRY_RAYTRACER_H_

#include <memory>
#include <vector>

#include "ion/math/vector.h"
#include "seurat/geometry/mesh.h"

namespace seurat {
namespace geometry {

// A wrapper around Embree for tracing single rays through static scenes of
// triangles.
class Raytracer {
 public:
  struct Intersection {
    // The t-value associated with the intersection, the amount to scale the ray
    // direction to reach the intersection point.
    float t_hit;

    // The index of the triangle, specified by the order of the triangle in the
    // index_buffer passed to "Build".
    int triangle_index;
  };

  // Custom destructor to free the embree scene.
  ~Raytracer();

  // Constructs a Raytracer to cast rays through a scene consisting of triangles
  // specified by the vertex and index buffers.
  static std::unique_ptr<Raytracer> Build(
      const std::vector<ion::math::Point3f>& vertex_buffer,
      const std::vector<int>& index_buffer);

  // Constructs a Raytracer to cast rays through the specified mesh.
  static std::unique_ptr<Raytracer> Build(const Mesh& mesh);

  // Casts a ray through the scene, returning whether any triangle was hit
  // within the specified |t_max| limit.
  //
  // If a triangle was hit, the t-value associated with the intersection is
  // returned via |t_hit|, and the index of the triangle, specified by the start
  // of the triangle in the index_buffer passed to "Build", is returned via
  // |triangle_index_hit|.
  bool FindFirstHit(const ion::math::Point3f& ray_origin,
                    const ion::math::Vector3f& ray_direction, float t_max,
                    float* t_hit, int* triangle_index_hit) const;

  // Casts an infinitely-long ray through the scene.
  bool FindFirstHit(const ion::math::Point3f& ray_origin,
                    const ion::math::Vector3f& ray_direction, float* t_hit,
                    int* triangle_index_hit) const {
    const float t_max = std::numeric_limits<float>::infinity();
    return FindFirstHit(ray_origin, ray_direction, t_max, t_hit,
                        triangle_index_hit);
  }

  // Casts a ray through the scene, returning all intersections.
  //
  // Note: Intersections may be returned in any order.
  void FindAllIntersections(const ion::math::Point3f& ray_origin,
                            const ion::math::Vector3f& ray_direction,
                            std::vector<Intersection>* intersections) const;

  // Counts the number of triangle intersections along the ray segment until
  // |t_max|.
  int CountIntersections(const ion::math::Point3f& ray_origin,
                         const ion::math::Vector3f& ray_direction, float t_max,
                         int max_intersection_count) const;

 private:
  // A struct containing an embree RTCScene.
  //
  // This PIMPL-style indirection is necessary to isolate problematic
  // preprocessor directives found in embree headers.
  struct RTCSceneHandle;

  explicit Raytracer(std::unique_ptr<RTCSceneHandle> scene);

  // A pointer to a struct containing the Embree scene built from triangles.
  std::unique_ptr<RTCSceneHandle> scene_;
};

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_RAYTRACER_H_
