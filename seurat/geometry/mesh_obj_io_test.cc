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

#include "seurat/geometry/mesh_obj_io.h"

#include <vector>

#include "gtest/gtest.h"
#include "seurat/geometry/mesh.h"

namespace seurat {
namespace geometry {

using ion::math::Point2f;
using ion::math::Point3f;
using Triangle = Mesh::Triangle;

TEST(MeshObjIo, WriteObj) {
  Mesh mesh(1);
  std::vector<Point3f> positions;
  positions.push_back(Point3f(0.0f, 0.0f, 0.0f));
  positions.push_back(Point3f(1.0f, 0.0f, 0.0f));
  positions.push_back(Point3f(1.0f, 1.0f, 0.0f));
  positions.push_back(Point3f(0.0f, 1.0f, 0.0f));
  std::vector<Point3f> tex_coords{
      Point3f(0.0f, 0.0f, 1.0f), Point3f(1.0f, 0.0f, 1.0f),
      Point3f(1.0f, 1.0f, 1.0f), Point3f(0.0f, 1.0f, 1.0f)};

  for (int i = 0; i < 4; ++i) {
    mesh.AppendVertex(positions[i], {{tex_coords[i]}});
  }
  std::vector<Triangle> triangles;
  triangles.push_back({{0, 1, 2}});
  triangles.push_back({{2, 3, 0}});
  for (const Triangle& triangle : triangles) {
    mesh.AppendTriangle(triangle);
  }

  std::string obj_data;
  MeshObjIo::WriteObj(&obj_data, mesh);
  std::string expected(
      "v 0.000000 0.000000 0.000000\n"  //
      "v 1.000000 0.000000 0.000000\n"  //
      "v 1.000000 1.000000 0.000000\n"  //
      "v 0.000000 1.000000 0.000000\n"  //
      "vt 0.000000 0.000000\n"          //
      "vt 1.000000 0.000000\n"          //
      "vt 1.000000 1.000000\n"          //
      "vt 0.000000 1.000000\n"          //
      "f 1/1 2/2 3/3\n"                 //
      "f 3/3 4/4 1/1\n");
  EXPECT_EQ(expected, obj_data);
}

}  // namespace geometry
}  // namespace seurat
