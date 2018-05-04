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

#include "seurat/baker/diffuse/diffuse_pipeline.h"

#include <algorithm>
#include <array>

#include "ion/gfx/image.h"
#include "ion/math/vector.h"
#include "ion/port/timer.h"
#include "seurat/artifact/artifact.h"
#include "seurat/baker/framework/frame.h"
#include "seurat/baker/framework/radiance_accumulator.h"
#include "seurat/baker/framework/ray_classifier.h"
#include "seurat/baker/framework/silhouette.h"
#include "seurat/baker/framework/silhouette_accumulator.h"
#include "seurat/base/color.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/base/parallel.h"
#include "seurat/base/util.h"
#include "seurat/image/filter.h"
#include "seurat/image/image_util.h"

namespace seurat {
namespace baker {
namespace diffuse {

using artifact::Artifact;
using image::Image1f;
using image::Image3f;
using image::Image4f;
using image::ImageView4f;
using ion::math::Matrix3f;
using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2f;
using ion::math::Vector2i;

base::Status DiffusePipeline::Process(Artifact* artifact) const {
  std::vector<Frame> frames;
  SEURAT_RETURN_IF_ERROR(frame_generator_->GenerateFrames(&frames));
  // Sort according to draw order.
  //
  // The order of the final quads is the same as the order of the frames.
  std::sort(frames.begin(), frames.end(),
            [](const Frame& lhs, const Frame& rhs) {
              return lhs.draw_order < rhs.draw_order;
            });
  const int num_frames = frames.size();

  std::vector<Vector2i> texture_size_per_frame(num_frames);
  texture_sizer_->ComputeTextureSizes(frames,
                                      absl::MakeSpan(texture_size_per_frame));

  std::vector<Image4f> textures(num_frames);
  for (int i = 0; i < num_frames; ++i) {
    textures[i].Resize(texture_size_per_frame[i]);
  }

  SEURAT_RETURN_IF_ERROR(
      baker_->BakeTextures(frames, absl::MakeSpan(textures)));

  base::ParallelFor(thread_count_, num_frames, [&](int frame_index) {
    color_processor_->ProcessColors(absl::MakeSpan(textures[frame_index]));
  });

  std::vector<geometry::IndexedQuad> rasterized_quads(num_frames);
  for (int i = 0; i < num_frames; ++i) {
    rasterized_quads[i] = {frames[i].quad, frames[i].texcoord_w, i};
  }

  *artifact = Artifact{};
  artifact->quad_mesh = std::make_shared<geometry::QuadMesh>(
      std::move(rasterized_quads), std::move(textures));
  return base::OkStatus();
}

}  // namespace diffuse
}  // namespace baker
}  // namespace seurat
