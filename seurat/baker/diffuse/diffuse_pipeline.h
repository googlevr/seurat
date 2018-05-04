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

#ifndef VR_SEURAT_BAKER_DIFFUSE_DIFFUSE_PIPELINE_H_
#define VR_SEURAT_BAKER_DIFFUSE_DIFFUSE_PIPELINE_H_

#include <memory>

#include "seurat/artifact/artifact.h"
#include "seurat/artifact/artifact_processor.h"
#include "seurat/baker/diffuse/diffuse_baker.h"
#include "seurat/baker/framework/frame_generator.h"
#include "seurat/baker/framework/texture_sizer.h"
#include "seurat/image/color_processor.h"
#include "seurat/mesh/mesh_component.h"

namespace seurat {
namespace baker {
namespace diffuse {

// Executes the given FrameGenerator & transforms the resulting Frames into a
// triangle mesh with an RGBA texture atlas.
class DiffusePipeline : public artifact::ArtifactProcessor {
 public:
  DiffusePipeline(int thread_count,
                  std::unique_ptr<FrameGenerator> frame_generator,
                  std::unique_ptr<TextureSizer> texture_sizer,
                  std::unique_ptr<DiffuseBaker> baker,
                  std::unique_ptr<image::ColorProcessor> color_processor)
      : thread_count_(thread_count),
        frame_generator_(std::move(frame_generator)),
        texture_sizer_(std::move(texture_sizer)),
        baker_(std::move(baker)),
        color_processor_(std::move(color_processor)) {}
  ~DiffusePipeline() override = default;

  // ArtifactProcessor implementation.
  base::Status Process(artifact::Artifact* artifact) const override;

 private:
  const int thread_count_;

  // Generates the collection of Frames used to model the geometry of the scene.
  const std::unique_ptr<FrameGenerator> frame_generator_;

  // Determines (initial) texture size.
  const std::unique_ptr<TextureSizer> texture_sizer_;

  // Rasterizes onto frames by performing a pass over all views.
  const std::unique_ptr<DiffuseBaker> baker_;

  // Post-processes color data, e.g. tone mapping.
  const std::unique_ptr<image::ColorProcessor> color_processor_;
};

}  // namespace diffuse
}  // namespace baker
}  // namespace seurat

#endif  // VR_SEURAT_BAKER_DIFFUSE_DIFFUSE_PIPELINE_H_
