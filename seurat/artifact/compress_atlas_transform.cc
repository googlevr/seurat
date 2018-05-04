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

#include "seurat/artifact/compress_atlas_transform.h"

#include "seurat/mesh/mesh_component_util.h"

namespace seurat {
namespace artifact {

base::Status CompressAtlasTransform::Process(Artifact* artifact) const {
  CHECK(artifact->texture)
      << "The input artifact does not support the required representation";
  // Extract and compress the atlas.
  const image::Image4f& texture_atlas = *artifact->texture;

  ion::gfx::ImagePtr ion_texture_atlas = codec_->Compress(texture_atlas);

  auto round_trip_atlas = std::make_shared<image::Image4f>();
  *round_trip_atlas = codec_->Decompress(ion_texture_atlas);

  artifact->texture = round_trip_atlas;
  artifact->component.reset(
      mesh::MeshComponentUtil::FromMesh(*artifact->mesh, ion_texture_atlas)
          .release());

  return base::OkStatus();
}

}  // namespace artifact
}  // namespace seurat
