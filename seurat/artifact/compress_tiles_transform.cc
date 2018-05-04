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

#include "seurat/artifact/compress_tiles_transform.h"

namespace seurat {
namespace artifact {

using geometry::QuadMesh;

base::Status CompressTilesTransform::Process(Artifact* artifact) const {
  CHECK(artifact->quad_mesh)
      << "The input artifact does not support the required representation";
  const QuadMesh& uncompressed = *artifact->quad_mesh;

  // Compress the tiles.
  std::vector<image::Image4f> compressed_textures(uncompressed.textures.size());
  codec_->Compress(uncompressed.textures, absl::MakeSpan(compressed_textures));

  artifact->quad_mesh = std::make_shared<QuadMesh>(
      artifact->quad_mesh->quads, std::move(compressed_textures));

  return base::OkStatus();
}

}  // namespace artifact
}  // namespace seurat
