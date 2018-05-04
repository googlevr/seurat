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

#ifndef VR_SEURAT_ARTIFACT_SEPARATE_OPAQUE_H_
#define VR_SEURAT_ARTIFACT_SEPARATE_OPAQUE_H_

#include <memory>
#include <utility>

#include "seurat/artifact/artifact.h"
#include "seurat/artifact/artifact_processor.h"
#include "seurat/base/status.h"

namespace seurat {
namespace artifact {

// Takes as input a QuadMesh artifact and output another QuadMesh artifact.
// The output QuadMesh consists of only those quads whose textures are
// either fully opaque (in case |retain_opaque| is set to true) or translucent
// (in case |retain_opaque| is set to false.
class SeparateOpaque : public ArtifactProcessor {
 public:
  enum class Retain { kRetainOpaque, kRetainTranslucent };
  SeparateOpaque(Retain retain, float alpha_threshold)
      : retain_(retain), alpha_threshold_(alpha_threshold) {}

  base::Status Process(Artifact* artifact) const override;

 private:
  // Defines which quads will be passed to the next artifact processing stage.
  Retain retain_;

  // Defines a threshold normalized to the interval [0.0f, 1.0f]. A quad is
  // deemed transparent if it has texels with an alpha below this threshold.
  float alpha_threshold_;
};

}  // namespace artifact
}  // namespace seurat

#endif  // VR_SEURAT_ARTIFACT_SEPARATE_OPAQUE_H_
