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

#ifndef VR_SEURAT_ARTIFACT_COMPRESS_TILES_TRANSFORM_H_
#define VR_SEURAT_ARTIFACT_COMPRESS_TILES_TRANSFORM_H_

#include <memory>
#include <utility>

#include "seurat/artifact/artifact.h"
#include "seurat/artifact/artifact_processor.h"
#include "seurat/base/status.h"
#include "seurat/compressor/rgba/rgba_compressor.h"
#include "seurat/image/atlaser.h"

namespace seurat {
namespace artifact {

// Compresses an Artifact's QuadMesh using a RgbaCompressor.
//
// Only the Artifact.quad_mesh is replaced.  All other elements are unmodified.
class CompressTilesTransform : public ArtifactProcessor {
 public:
  explicit CompressTilesTransform(
      std::unique_ptr<compressor::RgbaCompressor> codec,
      std::shared_ptr<image::Atlaser> atlaser)
      : codec_(std::move(codec)), atlaser_(std::move(atlaser)) {}

  // ArtifactProcessor implementation.
  base::Status Process(Artifact* artifact) const override;

 private:
  // Defines the compression operation for Transform.
  std::unique_ptr<compressor::RgbaCompressor> codec_;

  // Defines the process of laying out texture tiles in a texture atlas,
  // including possible constraints on the atlas size.
  std::shared_ptr<image::Atlaser> atlaser_;
};

}  // namespace artifact
}  // namespace seurat

#endif  // VR_SEURAT_ARTIFACT_COMPRESS_TILES_TRANSFORM_H_
