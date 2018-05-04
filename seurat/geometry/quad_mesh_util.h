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

#ifndef VR_SEURAT_GEOMETRY_QUAD_MESH_UTIL_H_
#define VR_SEURAT_GEOMETRY_QUAD_MESH_UTIL_H_

#include <memory>

#include "absl/types/span.h"
#include "seurat/geometry/mesh.h"
#include "seurat/geometry/quad.h"
#include "seurat/geometry/quad_mesh.h"
#include "seurat/image/atlaser.h"
#include "seurat/image/image.h"

namespace seurat {
namespace geometry {

// Atlasses a QuadMesh (set of textured quadrilaterals) into a Mesh & Texture.
//
// |mesh| or |atlas| may be nullptr if not required.
void MeshFromQuadMesh(absl::Span<const IndexedQuad> rasterized_quads,
                      absl::Span<const image::Image4f> textures,
                      const image::Atlaser& atlaser, geometry::Mesh* mesh,
                      image::Image4f* atlas);

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_QUAD_MESH_UTIL_H_
