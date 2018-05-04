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

#include "seurat/artifact/atlas_mesh_transform.h"

#include <memory>

#include "seurat/artifact/artifact.h"
#include "seurat/geometry/mesh.h"
#include "seurat/geometry/mesh_util.h"
#include "seurat/geometry/quad_mesh_util.h"
#include "seurat/image/image_util.h"
#include "seurat/mesh/mesh_component_util.h"

namespace seurat {
namespace artifact {

using geometry::Mesh;
using geometry::QuadMesh;
using image::Image4f;

base::Status AtlasMeshTransform::Process(Artifact* artifact) const {
  CHECK(artifact->quad_mesh);
  const QuadMesh& quad_mesh = *artifact->quad_mesh;

  constexpr int kTextureCount = 1;
  auto mesh = std::make_shared<Mesh>(kTextureCount);
  auto texture = std::make_shared<Image4f>();
  geometry::MeshFromQuadMesh(quad_mesh.quads, quad_mesh.textures, *atlaser_,
                             mesh.get(), texture.get());
  artifact->mesh = mesh;
  artifact->texture = texture;

  ion::gfx::ImagePtr ion_atlas =
      image::ConvertSeuratImageToIonImage(*artifact->texture);
  artifact->component.reset(
      mesh::MeshComponentUtil::FromMesh(*artifact->mesh, ion_atlas).release());

  return base::OkStatus();
}

}  // namespace artifact
}  // namespace seurat
