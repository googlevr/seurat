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

#include "seurat/artifact/mesh_artifact_util.h"

#include "seurat/base/color.h"
#include "seurat/geometry/mesh.h"
#include "seurat/geometry/plane.h"
#include "seurat/geometry/triangle.h"
#include "seurat/image/image.h"

namespace seurat {
namespace artifact {

using ion::math::Point3f;
using ion::math::Vector2i;

using base::Color4f;
using geometry::Mesh;
using geometry::Plane3f;
using geometry::Triangle3f;
using image::Image4f;

namespace {

Triangle3f Triangle3fFromMeshTriangle(const Mesh& mesh,
                                      const Mesh::Triangle& triangle) {
  const std::vector<ion::math::Point3f>& positions = mesh.GetPositions();
  return Triangle3f{
      {positions[triangle[0]], positions[triangle[1]], positions[triangle[2]]}};
}

void FlipTriangleToFaceOrigin(const Mesh& mesh, Mesh::Triangle* triangle) {
  Plane3f triangle_plane =
      geometry::PlaneFromTriangle(Triangle3fFromMeshTriangle(mesh, *triangle));
  // D < 0.0 means the origin and thus the headbox center is behind the plane of
  // the triangle, so flip it.
  if (triangle_plane.GetD() < 0.0f) {
    using std::swap;
    swap((*triangle)[0], (*triangle)[2]);
  }
}

}  // namespace

base::Status FlipMeshFacesTransform::Process(Artifact* artifact) const {
  CHECK(artifact->mesh);
  auto mutable_mesh = std::make_shared<Mesh>(*artifact->mesh);

  std::vector<Mesh::Triangle>& triangles = mutable_mesh->GetTriangles();
  for (Mesh::Triangle& triangle : triangles) {
    FlipTriangleToFaceOrigin(*mutable_mesh, &triangle);
  }

  artifact->mesh = mutable_mesh;
  return base::OkStatus();
}

}  // namespace artifact
}  // namespace seurat
