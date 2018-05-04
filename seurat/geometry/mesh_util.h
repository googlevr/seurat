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

#ifndef VR_SEURAT_GEOMETRY_MESH_UTIL_H_
#define VR_SEURAT_GEOMETRY_MESH_UTIL_H_

#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/base/util.h"
#include "seurat/geometry/mesh.h"

namespace seurat {
namespace geometry {

// Appends a triangle fan to the |mesh|.
//
// The added triangles tessellate the convex-polygon defined by |positions|
// and |tex_coords| such that each resulting triangle shares the vertex at
// positions.front().
//
// |tex_coords| must contain a texture-coordinate for each point in
// |positions| or be empty, resulting in zero-valued texture-coordinates for
// the triangles.
void AppendTriangleFan(
    absl::Span<const ion::math::Point3f> positions,
    absl::Span<const absl::Span<const ion::math::Point3f>> tex_coords,
    Mesh* mesh);

// Estimates the total area of all triangles of the given |mesh| if they were
// all projected onto a sphere centered at the specified point.
//
// The estimated area is returned as a multiple of the total area of the sphere
// (i.e. 1 indicates that the area of the mesh's projection is equal to the area
// of the sphere on which it is projected).
//
// The estimation assumes that the triangles are small enough that a small-angle
// approximation is applicable.  In other words, this assumes that the angle
// from a triangle vertex to |sphere_center| to another triangle vertex is small
// enough that sin(theta)=theta.
float EstimateProjectedArea(const ion::math::Point3f& sphere_center,
                            const Mesh& mesh);

// Transforms a Mesh into an equivalent Mesh for which the index-buffer is the
// identity (essentially unindexed).
Mesh ToUnindexedMesh(const Mesh& indexed);

// Splits the input mesh, as necessary, into meshes with the specified maximum
// number of vertices per mesh.
//
// This is useful for exporting meshes with only 16-bit index buffers.
//
// Triangle ordering is preserved in the output sequence of meshes.
std::vector<Mesh> SplitMeshByVertexCount(const Mesh& original,
                                         int max_vertex_count);

}  // namespace geometry
}  // namespace seurat

#endif  // VR_SEURAT_GEOMETRY_MESH_UTIL_H_
