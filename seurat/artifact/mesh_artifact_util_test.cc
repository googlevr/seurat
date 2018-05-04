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

#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/base/color.h"
#include "seurat/geometry/mesh.h"
#include "seurat/geometry/plane.h"
#include "seurat/geometry/triangle.h"
#include "seurat/image/image.h"

namespace seurat {
namespace artifact {
namespace {

using ion::math::Point3f;
using ion::math::Vector2i;

using base::Color4f;
using geometry::Mesh;
using geometry::Plane3f;
using geometry::Triangle3f;
using image::Image4f;

TEST(MeshArtifactUtilTest, FlipMeshFacesTransformProducesFrontFacingTriangles) {
  // 1. Builds the test mesh data: top-down view of triangles (abc) and (dae).
  //    Note that (abc) is backfacing as it has clockwise order viewed from the
  //    origin, while (dae) is counter-clockwise from the origin. This matches
  //    PlaneFromTriangle's definition of the positive halfspace of the plane.
  //    |
  //  d-+--ae
  //    |  |
  //----+--+------> x
  //    |  |
  //    |  bc
  //    |
  //    v +z
  const int kNoTextureCoordinates = 0;
  Mesh mesh(kNoTextureCoordinates);
  const Point3f kA(1.0f, 0.0f, -1.0f);
  const Point3f kB(1.0f, 1.0f, 1.0f);
  const Point3f kC(1.0f, 0.0f, 1.0f);
  const Point3f kD(-1.0f, 0.0f, -1.0f);
  const Point3f kE(1.0f, 1.0f, -1.0f);
  mesh.AppendVertex(kA, {});
  mesh.AppendVertex(kB, {});
  mesh.AppendVertex(kC, {});
  mesh.AppendVertex(kD, {});
  mesh.AppendVertex(kE, {});
  mesh.AppendTriangle(Mesh::Triangle({{0, 1, 2}}));
  mesh.AppendTriangle(Mesh::Triangle({{3, 0, 4}}));
  ASSERT_EQ(mesh.GetPositions().size(), 5);
  ASSERT_EQ(mesh.GetTriangles().size(), 2);

  // 2. Checks the input data contains both facing directions.
  const Point3f kOrigin = Point3f::Zero();
  bool found_forward = false;
  bool found_backward = false;
  for (Mesh::Triangle const& triangle : mesh.GetTriangles()) {
    const Triangle3f t = {{mesh.GetPositions()[triangle[0]],
                           mesh.GetPositions()[triangle[1]],
                           mesh.GetPositions()[triangle[2]]}};
    const Plane3f triangle_plane = geometry::PlaneFromTriangle(t);
    const float signed_distance_to_origin =
        triangle_plane.SignedDistanceToPoint(kOrigin);
    if (signed_distance_to_origin > 0.0f) {
      found_forward = true;
    }
    if (signed_distance_to_origin < 0.0f) {
      found_backward = true;
    }
  }
  ASSERT_TRUE(found_forward);
  ASSERT_TRUE(found_backward);

  // Builds the input artifact with inconsistent facing.
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  Artifact artifact;
  artifact.mesh = std::make_shared<Mesh>(mesh);
  artifact.texture = std::make_shared<Image4f>(Vector2i(1, 1), kRed);

  // 3. Produces a mesh artifact with triangles oriented toward the eye from the
  // input mesh artifact; the code under test.
  FlipMeshFacesTransform transform;
  EXPECT_TRUE(transform.Process(&artifact).ok());

  // 4. Verifies all triangles face the origin.
  const Mesh& output_mesh = *artifact.mesh;
  EXPECT_EQ(mesh.GetPositions().size(), output_mesh.GetPositions().size());
  EXPECT_EQ(mesh.GetTriangles().size(), output_mesh.GetTriangles().size());
  std::vector<Mesh::Triangle> output_triangles = output_mesh.GetTriangles();
  for (Mesh::Triangle triangle : output_triangles) {
    const Triangle3f t = {{output_mesh.GetPositions()[triangle[0]],
                           output_mesh.GetPositions()[triangle[1]],
                           output_mesh.GetPositions()[triangle[2]]}};
    const Plane3f triangle_plane = geometry::PlaneFromTriangle(t);
    const float signed_distance_to_origin =
        triangle_plane.SignedDistanceToPoint(kOrigin);
    EXPECT_GE(signed_distance_to_origin, 0.0f);
  }
}

}  // namespace
}  // namespace artifact
}  // namespace seurat
