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

#ifndef VR_SEURAT_ARTIFACT_ATLAS_MESH_TRANSFORM_H_
#define VR_SEURAT_ARTIFACT_ATLAS_MESH_TRANSFORM_H_

#include <memory>
#include <utility>

#include "seurat/artifact/artifact.h"
#include "seurat/artifact/artifact_processor.h"
#include "seurat/base/status.h"
#include "seurat/image/atlaser.h"

namespace seurat {
namespace artifact {

// Atlases a QuadMesh into an Image4f, Mesh, and Component.
class AtlasMeshTransform : public ArtifactProcessor {
 public:
  explicit AtlasMeshTransform(std::shared_ptr<image::Atlaser> atlaser)
      : atlaser_(std::move(atlaser)) {}

  base::Status Process(Artifact* artifact) const override;

 private:
  // Defines the process of laying out texture tiles in a texture atlas,
  // including possible constraints on the atlas size.
  std::shared_ptr<image::Atlaser> atlaser_;
};

}  // namespace artifact
}  // namespace seurat

#endif  // VR_SEURAT_ARTIFACT_ATLAS_MESH_TRANSFORM_H_
