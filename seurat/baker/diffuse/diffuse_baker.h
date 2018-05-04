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

#ifndef VR_SEURAT_BAKER_DIFFUSE_DIFFUSE_BAKER_H_
#define VR_SEURAT_BAKER_DIFFUSE_DIFFUSE_BAKER_H_

#include <memory>

#include "seurat/baker/framework/rasterizer.h"
#include "seurat/baker/framework/ray_classifier.h"
#include "seurat/image/filter.h"
#include "seurat/image/image.h"

namespace seurat {
namespace baker {
namespace diffuse {

// Bakes RGBA textures for a set of quads.
class DiffuseBaker {
 public:
  DiffuseBaker(int thread_count, std::unique_ptr<FrameRasterizer> rasterizer,
               float specular_filter_size,
               std::shared_ptr<image::Filter> pixel_filter)
      : thread_count_(thread_count),
        rasterizer_(std::move(rasterizer)),
        specular_filter_size_(specular_filter_size),
        pixel_filter_(std::move(pixel_filter)) {}

  // Bakes an RGBA texture for each piece of geometry.
  //
  // The resolution of the given |textures| is preserved and only pixel data is
  // modified.
  base::Status BakeTextures(absl::Span<const Frame> geometry,
                            absl::Span<image::Image4f> textures) const;

 private:
  // The number of threads to use.
  const int thread_count_;

  // Rasterizes onto frames by performing a pass over all views.
  const std::unique_ptr<FrameRasterizer> rasterizer_;

  // Determines how much specular highlights should be be baked in as if
  // sampling radiance from the origin.
  //
  // A smaller value will bake specular highlights which are visible from the
  // origin.
  //
  // A larger value yields a more diffuse-looking representation of the
  // scene.
  const float specular_filter_size_;

  // Filter used for pixel reconstruction.
  std::shared_ptr<image::Filter> pixel_filter_;
};

}  // namespace diffuse
}  // namespace baker
}  // namespace seurat

#endif  // VR_SEURAT_BAKER_DIFFUSE_DIFFUSE_BAKER_H_
