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

#include "seurat/baker/diffuse/diffuse_baker.h"

#include <algorithm>
#include <array>

#include "ion/gfx/image.h"
#include "ion/math/vector.h"
#include "seurat/baker/framework/frame.h"
#include "seurat/baker/framework/radiance_accumulator.h"
#include "seurat/baker/framework/ray_classifier.h"
#include "seurat/baker/framework/sample_accumulator.h"
#include "seurat/baker/framework/sample_accumulator_set.h"
#include "seurat/baker/framework/silhouette.h"
#include "seurat/baker/framework/silhouette_accumulator.h"
#include "seurat/base/color.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/base/parallel.h"
#include "seurat/base/util.h"
#include "seurat/geometry/mesh.h"
#include "seurat/image/filter.h"
#include "seurat/image/image_util.h"

namespace seurat {
namespace baker {
namespace diffuse {

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

base::Status DiffuseBaker::BakeTextures(
    absl::Span<const Frame> geometry,
    absl::Span<image::Image4f> textures) const {
  CHECK_EQ(geometry.size(), textures.size());
  const int num_frames = geometry.size();

  // Initialize accumulators to rasterize silhouettes & a radiance
  // representation.
  std::vector<std::shared_ptr<SilhouetteAccumulator>> silhouette_per_frame;
  std::vector<std::shared_ptr<RadianceAccumulator>> radiance_per_frame;
  for (int frame_index = 0; frame_index < num_frames; ++frame_index) {
    const Vector2i size = textures[frame_index].GetSize();
    // Use 8 * 8 = 64 samples per pixel when resolving silhouettes.
    const int kSilhouetteSuperSamplingFactor = 8;
    silhouette_per_frame.push_back(std::make_shared<SilhouetteAccumulator>(
        &geometry[frame_index], kSilhouetteSuperSamplingFactor,
        std::unique_ptr<SilhouetteBuffer>(new CompactSilhouetteBuffer(
            kSilhouetteSuperSamplingFactor * size))));

    // Resample radiance from the direction of the origin, assumed to be the
    // center of the viewing-volume.
    const Point3f eye = Point3f::Zero();

    // The variance of the gaussian filter to use when filtering a signal to be
    // sampled at 1-unit intervals.
    const float kSigmaUnit = 0.7f;
    // Standard deviation of the specular (directional) Gaussian filter.
    const float kSigmaEye = kSigmaUnit * specular_filter_size_;

    radiance_per_frame.push_back(std::make_shared<RadianceAccumulator>(
        &geometry[frame_index], eye, kSigmaEye, pixel_filter_, size));
  }

  // Wire up the accumulators into a single super-accumulator per Frame, so we
  // can feed this into the rasterizer.
  std::vector<std::shared_ptr<SampleAccumulator>> accumulators;
  for (int frame_index = 0; frame_index < num_frames; ++frame_index) {
    std::array<std::shared_ptr<SampleAccumulator>, 2> all_accumulators{
        {silhouette_per_frame[frame_index],  //
         radiance_per_frame[frame_index]}};
    accumulators.push_back(
        std::make_shared<SampleAccumulatorSet>(all_accumulators));
  }

  // Invoke the rasterizer, running a single pass over all RGBD Views.
  SEURAT_RETURN_IF_ERROR(rasterizer_->Run(geometry, accumulators));

  // Resolve all accumulators into a set of QuadMesh IndexedQuads.
  base::ParallelFor(thread_count_, num_frames, [&](int frame_index) {
    const Vector2i size = textures[frame_index].GetSize();

    Image1f silhouette(size);
    silhouette_per_frame[frame_index]->Resolve(&silhouette);

    Image3f rgb(size);
    radiance_per_frame[frame_index]->Resolve(&rgb);

    Image4f& rgba = textures[frame_index];
    for (int y = 0; y < size[1]; ++y) {
      for (int x = 0; x < size[0]; ++x) {
        rgba.At(x, y) = {rgb.At(x, y), silhouette.At(x, y)[0]};
      }
    }
  });

  return base::OkStatus();
}

}  // namespace diffuse
}  // namespace baker
}  // namespace seurat
