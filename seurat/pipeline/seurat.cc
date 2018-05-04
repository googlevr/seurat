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

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include "gflags/gflags.h"
#include "OpenEXR/IlmImf/ImfHeader.h"
#include "seurat/artifact/artifact.h"
#include "seurat/artifact/artifact_processor.h"
#include "seurat/base/progress.h"
#include "seurat/base/reporting.h"
#include "seurat/pipeline/pipeline.h"

using seurat::pipeline::Pipeline;

DEFINE_string(input_path, seurat::pipeline::proto::Flags().input_path(),
              "Path to the input manifest.json file.");
DEFINE_string(output_path, seurat::pipeline::proto::Flags().output_path(),
              "Base path to all output artifacts.");
DEFINE_string(cache_path, seurat::pipeline::proto::Flags().cache_path(),
              "Directory for all cached values.");
DEFINE_string(single_face, seurat::pipeline::proto::Flags().single_face(),
              "If not empty, process only the specified face. Must be one of "
              "'front', 'back', 'left', 'right', 'bottom', 'top'");
DEFINE_int32(triangle_count, seurat::pipeline::proto::Flags().triangle_count(),
             "The maximum number of triangles.");
DEFINE_double(overdraw_factor,
              seurat::pipeline::proto::Flags().overdraw_factor(),
              "The maximum amount of average overdraw.");
DEFINE_double(peak_overdraw_factor,
              seurat::pipeline::proto::Flags().peak_overdraw_factor(),
              "The maximum amount of overdraw from any view point.");
DEFINE_double(gamma, seurat::pipeline::proto::Flags().gamma(),
              "Gamma-correction exponent.");
DEFINE_double(specular_filter_size,
              seurat::pipeline::proto::Flags().specular_filter_size(),
              "The size of the filter used to bake specular highlights.  "
              "Smaller values bake sharper reflections.  Larger values blur "
              "these out, yielding a more diffuse-looking representation.");
DEFINE_bool(premultiply_alpha,
            seurat::pipeline::proto::Flags().premultiply_alpha(),
            "Determines whether output textures use premultiplied alpha.");
DEFINE_bool(evaluate, seurat::pipeline::proto::Flags().evaluate(),
            "Evalutes output");
DEFINE_double(ray_footprint, seurat::pipeline::proto::Flags().ray_footprint(),
              "The 'footprint' of a sample, along its depth.  Larger values "
              "help fill seams in the final geometry.");
DEFINE_double(pixels_per_degree,
              seurat::pipeline::proto::Flags().pixels_per_degree(),
              "Resolution of the target display in pixels per degree. This "
              "parameter is used to determine texture size.");
DEFINE_int32(texture_width, seurat::pipeline::proto::Flags().texture_width(),
             "The width of output textures.");
DEFINE_int32(texture_height, seurat::pipeline::proto::Flags().texture_height(),
             "The height of output textures.");
DEFINE_int32(texture_alignment,
             seurat::pipeline::proto::Flags().texture_alignment(),
             "Alignment constraint on individual texture tiles in the atlas.");
DEFINE_bool(content_adaptive_resolution,
            seurat::pipeline::proto::Flags().content_adaptive_resolution(),
            "Determines whether to adapt local texture resolution based on "
            "texture content.");
DEFINE_double(skybox_radius, seurat::pipeline::proto::Flags().skybox_radius(),
              "Half the side-length of the origin-centered skybox to clamp "
              "geometry. 0 indicates no skybox clamping should be performed.");
DEFINE_bool(fast_preview, seurat::pipeline::proto::Flags().fast_preview(),
            "Produces a fast, low-quality preview.");
DEFINE_bool(report_progress, true, "Print progress updates to stdout.");
DEFINE_bool(z_buffer, seurat::pipeline::proto::Flags().z_buffer(),
            "Whether the output will be rendered with z-test.");
DEFINE_bool(projective_texture_mapping, false,
            "Enables projective texture mapping. Otherwise, object space "
            "texture mapping is used.");
DEFINE_bool(
    separate_opaque, seurat::pipeline::proto::Flags().separate_opaque(),
    "Determines whether separate meshes and texture atlases will be output for"
    "opaque and translucent parts of the scene.");
DEFINE_double(alpha_threshold,
              seurat::pipeline::proto::Flags().alpha_threshold(),
              "Defines the threshold for deciding whether a texture is opaque "
              "or translucent. A value between 0.0 and 1.0 is expected.");
DEFINE_string(pixel_filter, seurat::pipeline::proto::Flags().pixel_filter(),
              "Pixel filter for texture generation. Must be one of 'box', "
              "'bspline', 'gaussian'.");
DEFINE_double(headbox_radius, seurat::pipeline::proto::Flags().headbox_radius(),
              "Radius of the desired headbox size.");

namespace seurat {
namespace pipeline {

using artifact::Artifact;
using artifact::ArtifactProcessor;

namespace {

int seurat_main(int argc, char** argv) {
  const bool kRemoveParsedFlags = true;
  gflags::ParseCommandLineFlags(&argc, &argv, kRemoveParsedFlags);

  // Prevent races during OpenEXR initialization.
  Imf::staticInitialize();

  base::GetProgress()->Enable(FLAGS_report_progress);

  proto::Flags flags;
  flags.set_input_path(FLAGS_input_path);
  flags.set_output_path(FLAGS_output_path);
  flags.set_cache_path(FLAGS_cache_path);

  flags.set_single_face(FLAGS_single_face);

  flags.set_triangle_count(FLAGS_triangle_count);
  flags.set_overdraw_factor(FLAGS_overdraw_factor);
  flags.set_peak_overdraw_factor(FLAGS_peak_overdraw_factor);
  flags.set_texture_width(FLAGS_texture_width);
  flags.set_texture_height(FLAGS_texture_height);
  flags.set_content_adaptive_resolution(FLAGS_content_adaptive_resolution);
  flags.set_texture_alignment(FLAGS_texture_alignment);

  flags.set_gamma(FLAGS_gamma);
  flags.set_specular_filter_size(FLAGS_specular_filter_size);
  flags.set_premultiply_alpha(FLAGS_premultiply_alpha);
  flags.set_ray_footprint(FLAGS_ray_footprint);
  flags.set_pixels_per_degree(FLAGS_pixels_per_degree);
  flags.set_skybox_radius(FLAGS_skybox_radius);
  flags.set_fast_preview(FLAGS_fast_preview);
  flags.set_z_buffer(FLAGS_z_buffer);
  flags.set_projective_texture_mapping(FLAGS_projective_texture_mapping);
  flags.set_separate_opaque(FLAGS_separate_opaque);
  flags.set_alpha_threshold(FLAGS_alpha_threshold);
  flags.set_pixel_filter(FLAGS_pixel_filter);
  flags.set_evaluate(FLAGS_evaluate);
  flags.set_headbox_radius(FLAGS_headbox_radius);

  Pipeline pipeline(flags);

  Pipeline::Factories factories = pipeline.AssemblePipelineFactories();

  std::shared_ptr<const ArtifactProcessor> processing_pipeline =
      factories.processing_pipeline();

  Artifact artifact;
  base::Status status = processing_pipeline->Process(&artifact);
  if (!status.ok()) {
    base::SeuratError(status.error_message());
    return EXIT_FAILURE;
  }

  base::SeuratInfo("Processing completed successfully. Exiting.");
  return EXIT_SUCCESS;
}

}  // namespace
}  // namespace pipeline
}  // namespace seurat

int main(int argc, char* argv[]) {
  return seurat::pipeline::seurat_main(argc, argv);
}
