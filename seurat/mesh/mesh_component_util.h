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

#ifndef VR_SEURAT_MESH_MESH_COMPONENT_UTIL_H_
#define VR_SEURAT_MESH_MESH_COMPONENT_UTIL_H_

#include <memory>

#include "ion/gfx/image.h"
#include "ion/math/vector.h"
#include "seurat/geometry/mesh.h"
#include "seurat/mesh/mesh_component.h"

namespace seurat {
namespace mesh {

class MeshComponentUtil {
 public:
  // Sorting order for SortTrianglesByDistanceToEye method.
  enum class SortingOrder { kBackToFront = 0, kFrontToBack };

  // Constructs a unit-sized, axis-aligned, origin-centered cube, textured with
  // each face displaying the specified image.
  //
  // The exterior of the cube has faces oriented counter-clockwise.
  static std::unique_ptr<const MeshComponent> CreateCube(
      ion::gfx::ImagePtr cube_face_texture);

  // Converts a MeshComponent into a Mesh.
  static geometry::Mesh ToMesh(const MeshComponent& component);

  // Converts a Mesh into a MeshComponent using the specified
  // |texture_atlas|.
  static std::unique_ptr<const MeshComponent> FromMesh(
      const geometry::Mesh& mesh, ion::gfx::ImagePtr texture_atlas);

  // Sorts triangles in the |original| mesh by the distance from their center to
  // the given |eye| point.  The returned mesh will have triangles sorted as
  // specified by |sorting_order|.
  static std::unique_ptr<const MeshComponent> SortTrianglesByDistanceToEye(
      const MeshComponent& original, const ion::math::Point3f& eye,
      SortingOrder sorting_order);

  // Counts the edges between transparent & opaque texels.
  // This may be used to estimate mesh quality.
  //
  // Edges along the perimeter of the atlas are *not* counted.
  //
  // For example, we should expect that a billboard-cloud mesh consisting of
  // textures with sparse points (resulting in more edges between transparent &
  // opaque texels), will be less efficient to compress & render than one with
  // dense/compact textures.
  static int CountTransparentTexelEdges(const MeshComponent& mesh);
};

}  // namespace mesh
}  // namespace seurat

#endif  // VR_SEURAT_MESH_MESH_COMPONENT_UTIL_H_
