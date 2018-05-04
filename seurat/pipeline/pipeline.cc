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

#include "seurat/pipeline/pipeline.h"

#include <cstdlib>
#include <memory>
#include <string>

#include "absl/strings/str_cat.h"
#include "seurat/artifact/artifact.h"
#include "seurat/artifact/artifact_processor.h"
#include "seurat/artifact/atlas_mesh_transform.h"
#include "seurat/artifact/compress_atlas_transform.h"
#include "seurat/artifact/compress_tiles_transform.h"
#include "seurat/artifact/evaluation/artifact_evaluator.h"
#include "seurat/artifact/evaluation/cost_estimator.h"
#include "seurat/artifact/exr_exporter.h"
#include "seurat/artifact/ice_exporter.h"
#include "seurat/artifact/mesh_artifact_util.h"
#include "seurat/artifact/obj_exporter.h"
#include "seurat/artifact/png_exporter.h"
#include "seurat/artifact/separate_opaque.h"
#include "seurat/artifact/sort_atlas_tiles_transform.h"
#include "seurat/baker/diffuse/diffuse_baker.h"
#include "seurat/baker/diffuse/diffuse_pipeline.h"
#include "seurat/baker/framework/frame_generator.h"
#include "seurat/baker/framework/frame_sorter.h"
#include "seurat/baker/framework/rasterizer.h"
#include "seurat/baker/framework/ray_classifier.h"
#include "seurat/baker/framework/texture_sizer.h"
#include "seurat/base/parallel.h"
#include "seurat/base/reporting.h"
#include "seurat/component/component.h"
#include "seurat/compressor/rgba/rgba_compressor.h"
#include "seurat/compressor/rgba/rgba_compressor_util.h"
#include "seurat/compressor/rgba/rgba_nop_compressor.h"
#include "seurat/geometry/binning_point_cloud_builder.h"
#include "seurat/geometry/cube_face.h"
#include "seurat/geometry/point_cloud_builder.h"
#include "seurat/image/atlaser.h"
#include "seurat/image/color_processor.h"
#include "seurat/image/filter.h"
#include "seurat/image/fixed_width_atlaser.h"
#include "seurat/image/image.h"
#include "seurat/image/rgbaui8_codec.h"
#include "seurat/ingest/pattern_loader.h"
#include "seurat/ingest/point_cloud_assembler.h"
#include "seurat/ingest/point_cloud_assembler_factory.h"
#include "seurat/ingest/tone_mapping_view_group_loader.h"
#include "seurat/ingest/view_group_loader.h"
#include "seurat/ingest/view_group_loader_factory.h"
#include "seurat/tiler/tiler.h"

namespace seurat {
namespace pipeline {

using artifact::ArtifactEvaluator;
using artifact::ArtifactProcessor;
using artifact::ArtifactProcessorGroup;
using artifact::ArtifactProcessorSequence;
using artifact::EvaluationExporter;
using artifact::GeometryDistortionEvaluator;
using artifact::OverdrawEvaluator;
using artifact::RenderEvaluator;
using artifact::RenderSim;
using artifact::SeparateOpaque;
using baker::DiskCachingFrameGenerator;
using baker::FrameGenerator;
using baker::diffuse::DiffusePipeline;
using base::FileSystem;
using component::Component;
using compressor::RgbaCompressor;
using compressor::RgbaNopCompressor;
using geometry::CubeFace;
using image::Atlaser;
using image::Codec;
using image::ColorProcessor;
using image::ColorProcessorPipeline;
using image::FixedWidthAtlaser;
using image::GammaToneMapper;
using image::PremultipliedAlphaConverter;
using ingest::PatternLoader;
using ingest::PointCloudAssembler;
using ingest::PointCloudAssemblerFactory;
using ingest::ViewGroupLoader;
using ingest::ViewGroupLoaderFactory;
using ion::math::Vector2i;
using ion::math::Vector3i;
using tiler::TilerFactory;

namespace {

int GetThreadCount() { return base::GetNumberOfHardwareThreads(); }

}  // namespace

Pipeline::Factories Pipeline::AssemblePipelineFactories() const {
  // This function wires together the end-to-end Seurat pipeline.
  //
  // The returned Factories are entry points into various trees of functions
  // which define how to lazily-construct whatever service is to be created
  // (i.e. an ArtifactProcessingPipeline).
  //
  // Note that multiple redundant versions of the same dependencies might be
  // constructed (e.g. multiple identical ViewGroupLoader objects may be
  // constructed as dependencies of different systems), but this is usually
  // okay.

  // Validate input arguments.
  if (flags_.output_path().empty()) {
    base::SeuratFatal("Output path required");
  }
  if (flags_.input_path().empty() && flags_.pattern_name().empty()) {
    base::SeuratFatal("Input path required");
  }
  if (!flags_.single_face().empty() &&
      !geometry::IsCubeFaceString(flags_.single_face())) {
    base::SeuratFatal(
        "Face name must be empty or one of "  //
        "'front', 'back', 'left', 'right', 'bottom', 'top'");
  }

  const auto input_path_split = FileSystem::SplitPath(flags_.input_path());
  auto input_fs = std::make_shared<FileSystem>(input_path_split.first);
  const std::string input_file(input_path_split.second);

  const auto output_path_split = FileSystem::SplitPath(flags_.output_path());
  auto output_fs = std::make_shared<FileSystem>(output_path_split.first);
  const std::string output_file(output_path_split.second);

  const bool fast_geometry = flags_.fast_preview();
  // 8 views for fast preview.
  const int max_views =
      flags_.fast_preview() ? 8 : std::numeric_limits<int>::max();

  auto build_tiler = [=]() {
    TilerFactory::Parameters parameters;
    parameters.tile_count = flags_.triangle_count() / 2;
    parameters.thread_count = base::GetNumberOfHardwareThreads();
    parameters.min_subdivision_level = 3;
    parameters.max_subdivision_level = 7;
    parameters.overdraw_factor = flags_.overdraw_factor();
    parameters.peak_overdraw_factor = flags_.peak_overdraw_factor();
    // Use 2 less than 128, so that the total constraint vector (including
    // constraints for triangle-count and average overdraw) is cache aligned.
    const int kPeakOverdrawSamples = 126;
    parameters.peak_overdraw_samples = kPeakOverdrawSamples;
    parameters.headbox_radius = flags_.headbox_radius();
    parameters.skybox_radius = flags_.skybox_radius();
    parameters.fast = fast_geometry;
    return TilerFactory::CreateDefaultTiler(parameters);
  };

  PatternLoader::Parameters pattern_parameters{
      flags_.pattern_name(), flags_.pattern_image_size(),
      flags_.pattern_feature_size(), static_cast<float>(flags_.pattern_angle()),
      static_cast<float>(flags_.pattern_depth())};

  auto build_view_loader = [=]() {
    ViewGroupLoaderFactory::Parameters parameters;
    parameters.thread_count = GetThreadCount();
    parameters.single_face = flags_.single_face();
    parameters.max_views = max_views;
    parameters.skybox_radius = flags_.skybox_radius();
    parameters.pattern_parameters = pattern_parameters;
    std::unique_ptr<ViewGroupLoader> view_group_loader;
    ViewGroupLoaderFactory factory(parameters, input_fs);
    base::Status status =
        factory.CreateViewGroupLoader(input_file, &view_group_loader);
    if (!status.ok()) {
      base::SeuratFatal(status.error_message());
    }
    return view_group_loader;
  };

  auto build_point_cloud_assembler = [=]() {
    PointCloudAssemblerFactory::Parameters parameters;
    parameters.thread_count = GetThreadCount();
    parameters.pixels_per_degree = flags_.pixels_per_degree();
    return PointCloudAssemblerFactory(parameters)
        .CreatePointCloudAssembler(build_view_loader());
  };

  auto build_frame_generator = [=]() {
    std::unique_ptr<FrameGenerator> generator;
    generator.reset(new baker::TilerFrameGenerator(
        build_point_cloud_assembler(), build_tiler(),
        std::unique_ptr<baker::FrameSorter>(
            new baker::FrameSorter(GetThreadCount())),
        flags_.projective_texture_mapping()
            ? baker::TilerFrameGenerator::TextureParameterization::kProjective
            : baker::TilerFrameGenerator::TextureParameterization::
                  kObjectSpace));
    if (!flags_.cache_path().empty()) {
      // Wrap with a DiskCachingFrameGenerator.
      const char kCacheFilename[] = "geometry_cache.bin";
      std::shared_ptr<FileSystem> cache_fs;
      cache_fs.reset(new FileSystem(flags_.cache_path()));
      generator.reset(new DiskCachingFrameGenerator(kCacheFilename, cache_fs,
                                                    std::move(generator)));
    }
    return generator;
  };

  auto build_color_processor = [=]() {
    // Apply gamma-correction & then convert to premultiplied alpha (if
    // requested).
    std::vector<std::unique_ptr<ColorProcessor>> processors;
    processors.push_back(
        std::unique_ptr<ColorProcessor>(new GammaToneMapper(flags_.gamma())));
    if (flags_.premultiply_alpha()) {
      processors.push_back(
          std::unique_ptr<ColorProcessor>(new PremultipliedAlphaConverter()));
    }
    return std::unique_ptr<ColorProcessor>(
        new ColorProcessorPipeline(std::move(processors)));
  };

  auto build_atlaser = [=]() {
    return std::unique_ptr<image::Atlaser>(new image::FixedWidthAtlaser(
        Vector2i(flags_.texture_width(), flags_.texture_height())));
  };

  auto build_atlas_codec = [=]() {

    return std::unique_ptr<Codec>(new image::RgbaUi8Codec);
  };

  auto build_texture_sizer = [=]() {
    Vector2i block_size(flags_.texture_alignment(), flags_.texture_alignment());
    // Constructs a TextureSizer using the given scale factor.
    auto make_sizer_for_scale = [=](float scale) {
      std::unique_ptr<baker::TextureSizer> projected_area_texture_sizer(
          new baker::ProjectedAreaTextureSizer(scale *
                                               flags_.pixels_per_degree()));
      return std::unique_ptr<baker::TextureSizer>(new baker::BucketTextureSizer(
          std::move(projected_area_texture_sizer), block_size));
    };
    if (!flags_.content_adaptive_resolution()) {
      // If not bounding texture resolution with rate-distortion optimization,
      // limit it here in the TextureSizer by rescaling appropriately.
      return std::unique_ptr<baker::TextureSizer>(
          new baker::ConstrainedAtlasTextureSizer(build_atlaser(),
                                                  make_sizer_for_scale));
    } else {
      return make_sizer_for_scale(1.0f);
    }
  };

  auto build_ray_classifier = [=]() {
    // Rasterize with overlapping samples.
    //
    // Larger values result in greater duplication of texture information
    // among nearby frames. Smaller values may result in cracks due to
    // insufficient inpainting.
    const float kSecondaryFrameThreshold = flags_.ray_footprint();

    // The output artifacts are designed to be rendered with alpha-blending and
    // no z-test.
    baker::ProjectingRayClassifier::RenderingMode rendering_mode =
        baker::ProjectingRayClassifier::RenderingMode::kDrawOrder;
    if (flags_.z_buffer()) {
      rendering_mode = baker::ProjectingRayClassifier::RenderingMode::kZBuffer;
    }

    std::unique_ptr<baker::RayClassifier> classifier(
        new baker::ProjectingRayClassifier(GetThreadCount(), rendering_mode,
                                           kSecondaryFrameThreshold));

    return classifier;
  };

  auto build_pixel_filter = [=]() -> std::shared_ptr<image::Filter> {
    std::shared_ptr<image::Filter> filter;
    if (flags_.pixel_filter() == "box") {
      return std::make_shared<image::BoxFilter>();
    } else if (flags_.pixel_filter() == "bspline") {
      return std::make_shared<image::BSplineFilter>();
    } else if (flags_.pixel_filter() == "gaussian") {
      const float kSigma = 0.3f;
      const float kRadius = 1.5f;
      return std::make_shared<image::GaussianFilter>(kSigma, kRadius);
    } else {
      base::SeuratFatal(
          "pixel_filter must be one of: 'box', 'bspline', 'gaussian'.");
    }
  };

  auto build_rasterizer = [=]() {
    return std::unique_ptr<baker::FrameRasterizer>(new baker::FrameRasterizer(
        GetThreadCount(), build_view_loader(), build_ray_classifier()));
  };

  auto build_texture_compressor = [=]() {
    const Vector2i compressor_alignment(flags_.texture_alignment(),
                                        flags_.texture_alignment());
    return compressor::BuildAtlasSizeTargetRgbaCompressor(
        GetThreadCount(), compressor_alignment, build_atlaser());
  };

  // The following defines Artifact transformations to extend the pipeline.
  auto build_compress_tiles = [=]() {
    return std::make_shared<const artifact::CompressTilesTransform>(
        build_texture_compressor(), build_atlaser());
  };
  auto build_sort_atlas_tiles = [=]() {
    return std::make_shared<const artifact::SortAtlasTilesTransform>();
  };
  auto build_atlas_mesh = [=]() {
    return std::make_shared<const artifact::AtlasMeshTransform>(
        build_atlaser());
  };
  auto build_flip_mesh_faces = [=]() {
    return std::make_shared<const artifact::FlipMeshFacesTransform>();
  };
  auto build_compress_atlas = [=]() {
    return std::make_shared<const artifact::CompressAtlasTransform>(
        build_atlas_codec());
  };

  // This builds the stage that generates the initial artifact.
  auto build_diffuse_source = [=]() {
    // N.B. uses shared_ptr constructor instead of make_shared for return type
    // deduction for the lambda versus the other pipeline lambdas.
    return std::shared_ptr<const ArtifactProcessor>(new DiffusePipeline(
        GetThreadCount(), build_frame_generator(), build_texture_sizer(),
        std::unique_ptr<baker::diffuse::DiffuseBaker>(
            new baker::diffuse::DiffuseBaker(
                GetThreadCount(), build_rasterizer(),
                flags_.specular_filter_size(), build_pixel_filter())),
        build_color_processor()));
  };

  auto build_source_stage = [=]() {
    if (flags_.baker_type() == "diffuse") {
      return build_diffuse_source();
    } else {
      base::SeuratFatal("Unknown baker type: " + flags_.baker_type());
    }
  };

  auto build_evaluator = [=]() {
    auto build_tone_mapped_view_loader = [=]() {
      return std::unique_ptr<ViewGroupLoader>(
          new ingest::ToneMappingViewGroupLoader(build_color_processor(),
                                                 build_view_loader()));
    };

    std::vector<std::unique_ptr<ArtifactEvaluator>> evaluators;
    evaluators.emplace_back(
        new OverdrawEvaluator(build_view_loader(), GetThreadCount()));
    evaluators.emplace_back(
        new RenderEvaluator(GetThreadCount(),
                            std::unique_ptr<RenderSim>(new RenderSim(
                                GetThreadCount(), flags_.z_buffer())),
                            build_tone_mapped_view_loader()));
    evaluators.emplace_back(
        new GeometryDistortionEvaluator(GetThreadCount(), build_view_loader()));
    return std::shared_ptr<ArtifactProcessor>(
        new EvaluationExporter(output_fs, output_file, std::move(evaluators)));
  };

  auto build_exporter = [=](std::string suffix) {
    std::vector<std::shared_ptr<const ArtifactProcessor>> stages;
    stages.emplace_back(build_sort_atlas_tiles());
    stages.emplace_back(build_atlas_mesh());
    stages.emplace_back(build_flip_mesh_faces());
    stages.emplace_back(build_compress_atlas());
    std::vector<std::shared_ptr<const ArtifactProcessor>> exporters;
    std::string full_output_file = absl::StrCat(output_file, suffix);
    exporters.emplace_back(
        new artifact::IceExporter(output_fs, full_output_file));
    exporters.emplace_back(
        new artifact::ObjExporter(output_fs, full_output_file));
    exporters.emplace_back(
        new artifact::PngExporter(output_fs, full_output_file));
    exporters.emplace_back(
        new artifact::ExrExporter(output_fs, full_output_file));
    if (flags_.evaluate()) {
      exporters.push_back(build_evaluator());
    }
    stages.emplace_back(
        std::make_shared<const ArtifactProcessorGroup>((exporters)));

    return std::make_shared<const ArtifactProcessorSequence>(stages);
  };

  const std::string kNoSuffix = "";

  auto build_separate_opaque_pipeline = [=]() {
    std::array<std::pair<SeparateOpaque::Retain, std::string>, 2>
        parameter_sets = {
            {{SeparateOpaque::Retain::kRetainOpaque, "_opaque"},
             {SeparateOpaque::Retain::kRetainTranslucent, "_translucent"}}};
    std::vector<std::shared_ptr<const ArtifactProcessor>> branches;
    for (const auto& parameter_set : parameter_sets) {
      std::vector<std::shared_ptr<const ArtifactProcessor>> branch_stages;
      branch_stages.emplace_back(std::shared_ptr<const ArtifactProcessor>(
          new SeparateOpaque(parameter_set.first,
                             static_cast<float>(flags_.alpha_threshold()))));
      branch_stages.emplace_back(build_exporter(parameter_set.second));
      branches.emplace_back(
          std::make_shared<const ArtifactProcessorSequence>(branch_stages));
    }

    // Add a branch that will also output the complete mesh.
    branches.emplace_back(build_exporter(kNoSuffix));

    return std::make_shared<const ArtifactProcessorGroup>(branches);
  };

  auto build_processing_pipeline = [=]() {
    // Defines the sequence of pipeline stages.
    std::vector<std::shared_ptr<const ArtifactProcessor>> stages;
    stages.push_back(build_source_stage());
    if (flags_.content_adaptive_resolution()) {
      stages.push_back(build_compress_tiles());
    }
    if (flags_.separate_opaque()) {
      stages.push_back(build_separate_opaque_pipeline());
    } else {
      stages.push_back(build_exporter(kNoSuffix));
    }

    return std::shared_ptr<const ArtifactProcessor>(
        new ArtifactProcessorSequence(stages));
  };

  return {build_processing_pipeline};
}

}  // namespace pipeline
}  // namespace seurat
