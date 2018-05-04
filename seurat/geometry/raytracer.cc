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

#include "seurat/geometry/raytracer.h"

#include <algorithm>
#include <iterator>
#include <limits>

#include "include/embree2/rtcore.h"
#include "include/embree2/rtcore_ray.h"
#include "seurat/base/util.h"

namespace seurat {
namespace geometry {

namespace {

// Custom extended ray struct for use with Embree.
//
// This follows the pattern outlined in the Embree tutorials for adding custom
// data to an RTCRay by declaring a custom structure with the exact layout, but
// with additional members appended.
struct RTCORE_ALIGN(16) RTCRay2 {
  // Fields part of Embree's RTCRay.  See rtcore_ray.h.
  float org[3];
  float align0;

  float dir[3];
  float align1;

  float tnear;
  float tfar;

  float time;
  int mask;

  float ng[3];
  float align2;

  float u;
  float v;

  int geom_id;
  int prim_id;
  int inst_id;

  // Custom fields:

  // The number of intersected primitives.
  int intersect_count;

  // The maximum number of intersections to count before stopping traversal of
  // the ray.
  int max_intersect_count;

  // A pointer to a vector to store all (t hit, triangle index) pairs, or
  // nullptr if they should be ignored.
  std::vector<Raytracer::Intersection>* intersections;
};

// This must match the embree signature for an embree
// IntersectionFilterFunction.
void IntersectionFilter(void* ptr,
                        RTCRay2& ray) {  // NOLINT(runtime/references)
  ray.intersect_count++;
  if (ray.intersect_count < ray.max_intersect_count) {
    // Continue traversal.
    ray.geom_id = RTC_INVALID_GEOMETRY_ID;
  }
  if (ray.intersections) ray.intersections->push_back({ray.tfar, ray.prim_id});
}

}  // namespace

using ion::math::Point3f;
using ion::math::Vector3f;

struct Raytracer::RTCSceneHandle {
  RTCScene scene;
};

Raytracer::Raytracer(std::unique_ptr<RTCSceneHandle> scene)
    : scene_(std::move(scene)) {}

Raytracer::~Raytracer() { rtcDeleteScene(scene_->scene); }

std::unique_ptr<Raytracer> Raytracer::Build(
    const std::vector<Point3f>& vertex_buffer,
    const std::vector<int>& index_buffer) {
  CHECK_EQ(0, index_buffer.size() % 3);

  // Only use a single RTCDevice for all instances of the Raytracer class.
  static RTCDevice s_device = rtcNewDevice(nullptr);

  RTCScene scene = rtcDeviceNewScene(
      s_device, RTC_SCENE_STATIC | RTC_SCENE_HIGH_QUALITY, RTC_INTERSECT1);

  int num_vertices = vertex_buffer.size();
  int num_triangles = index_buffer.size() / 3;

  unsigned int mesh = rtcNewTriangleMesh(scene, RTC_GEOMETRY_STATIC,
                                         num_triangles, num_vertices);

  struct EmbreeVertex {
    float x, y, z;
    // Embree uses structs of 4 floats for each vertex for alignment.
    //
    // So, this last value is for padding only.
    float r;
  };
  EmbreeVertex* vertices =
      static_cast<EmbreeVertex*>(rtcMapBuffer(scene, mesh, RTC_VERTEX_BUFFER));
  for (int i = 0; i < num_vertices; ++i) {
    vertices[i].x = vertex_buffer[i][0];
    vertices[i].y = vertex_buffer[i][1];
    vertices[i].z = vertex_buffer[i][2];
    vertices[i].r = 1.0f;
  }
  rtcUnmapBuffer(scene, mesh, RTC_VERTEX_BUFFER);

  int* triangles =
      static_cast<int*>(rtcMapBuffer(scene, mesh, RTC_INDEX_BUFFER));
  std::copy(index_buffer.begin(), index_buffer.end(), triangles);
  rtcUnmapBuffer(scene, mesh, RTC_INDEX_BUFFER);

  // Set our custom intersection filter for rtcOccluded() queries.
  // This is used for CountIntersections()
  //
  // TODO(puneetl):  Allow disabling this when we do not need to use
  // CountIntersections() and wish to use rtcOccluded() for other uses.
  rtcSetOcclusionFilterFunction(
      scene, mesh, reinterpret_cast<RTCFilterFunc>(&IntersectionFilter));

  // Commit (build BVH) in single threaded mode to avoid issues with Embree's
  // internal tasking system on Windows.
  //
  // TODO(b/34642024): Switch back to multi-threaded when Embree was upgraded to
  // 2.13.0 and TBB 2017.
  rtcCommitThread(scene, 0, 1);

  return std::unique_ptr<Raytracer>(new Raytracer(
      std::unique_ptr<RTCSceneHandle>(new RTCSceneHandle{scene})));
}

std::unique_ptr<Raytracer> Raytracer::Build(const Mesh& mesh) {
  std::vector<Point3f> vertex_buffer;
  std::vector<int> index_buffer;
  for (const auto& position : mesh.GetPositions()) {
    vertex_buffer.push_back(position);
  }
  for (const auto& triangle : mesh.GetTriangles()) {
    std::copy(triangle.begin(), triangle.end(),
              std::back_inserter(index_buffer));
  }
  return Build(vertex_buffer, index_buffer);
}

bool Raytracer::FindFirstHit(const Point3f& ray_origin,
                             const Vector3f& ray_direction, float t_max,
                             float* t_hit, int* triangle_index_hit) const {
  RTCRay ray;
  ray.org[0] = ray_origin[0];
  ray.org[1] = ray_origin[1];
  ray.org[2] = ray_origin[2];
  ray.align0 = 1.0f;
  ray.dir[0] = ray_direction[0];
  ray.dir[1] = ray_direction[1];
  ray.dir[2] = ray_direction[2];
  ray.align1 = 1.0f;
  ray.tnear = 0.0f;
  ray.tfar = t_max;
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.instID = RTC_INVALID_GEOMETRY_ID;
  ray.mask = 0xFFFFFFFF;
  ray.time = 0.0f;

  rtcIntersect(scene_->scene, ray);

  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) {
    return false;
  }

  *t_hit = ray.tfar;
  *triangle_index_hit = ray.primID;
  return true;
}

void Raytracer::FindAllIntersections(
    const ion::math::Point3f& ray_origin,
    const ion::math::Vector3f& ray_direction,
    std::vector<Intersection>* intersections) const {
  intersections->clear();
  RTCRay2 ray;
  ray.org[0] = ray_origin[0];
  ray.org[1] = ray_origin[1];
  ray.org[2] = ray_origin[2];
  ray.align0 = 1.0f;
  ray.dir[0] = ray_direction[0];
  ray.dir[1] = ray_direction[1];
  ray.dir[2] = ray_direction[2];
  ray.align1 = 1.0f;
  ray.tnear = 0.0f;
  ray.tfar = std::numeric_limits<float>::infinity();  // No limit.
  ray.geom_id = RTC_INVALID_GEOMETRY_ID;
  ray.prim_id = RTC_INVALID_GEOMETRY_ID;
  ray.inst_id = RTC_INVALID_GEOMETRY_ID;
  ray.mask = 0xFFFFFFFF;
  ray.time = 0.0f;

  ray.intersect_count = 0;                                    // Ignored.
  ray.max_intersect_count = std::numeric_limits<int>::max();  // No limit.
  ray.intersections = intersections;  // Log intersections.

  rtcOccluded(scene_->scene, reinterpret_cast<RTCRay&>(ray));
}

int Raytracer::CountIntersections(const ion::math::Point3f& ray_origin,
                                  const ion::math::Vector3f& ray_direction,
                                  float t_max,
                                  int max_intersection_count) const {
  RTCRay2 ray;
  ray.org[0] = ray_origin[0];
  ray.org[1] = ray_origin[1];
  ray.org[2] = ray_origin[2];
  ray.align0 = 1.0f;
  ray.dir[0] = ray_direction[0];
  ray.dir[1] = ray_direction[1];
  ray.dir[2] = ray_direction[2];
  ray.align1 = 1.0f;
  ray.tnear = 0.0f;
  ray.tfar = t_max;  // Limit t_max.
  ray.geom_id = RTC_INVALID_GEOMETRY_ID;
  ray.prim_id = RTC_INVALID_GEOMETRY_ID;
  ray.inst_id = RTC_INVALID_GEOMETRY_ID;
  ray.mask = 0xFFFFFFFF;
  ray.time = 0.0f;

  ray.intersect_count = 0;
  ray.max_intersect_count =
      max_intersection_count;   // Limit intersection count.
  ray.intersections = nullptr;  // Don't log intersections.

  rtcOccluded(scene_->scene, reinterpret_cast<RTCRay&>(ray));

  return ray.intersect_count;
}

}  // namespace geometry
}  // namespace seurat
