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

#ifndef VR_SEURAT_BAKER_FRAMEWORK_RADIANCE_ACCUMULATOR_H_
#define VR_SEURAT_BAKER_FRAMEWORK_RADIANCE_ACCUMULATOR_H_

#include <memory>

#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/baker/framework/frame.h"
#include "seurat/baker/framework/ray_bundle.h"
#include "seurat/baker/framework/sample_accumulator.h"
#include "seurat/base/color.h"
#include "seurat/image/filter.h"
#include "seurat/image/image.h"

namespace seurat {
namespace baker {

// A SampleAccumulator which accumulates Radiance towards a specified eye point.
//
// Samples are filtered by their associated rays' distance to the eye-point,
// using a gaussian filter with the specified sigma value, sigma_eye.
class RadianceAccumulator : public SampleAccumulator {
 public:
  RadianceAccumulator(const Frame* frame, const ion::math::Point3f& eye,
                      float sigma_eye,
                      std::shared_ptr<image::Filter> pixel_filter,
                      const ion::math::Vector2i& resolution)
      : frame_(frame),
        eye_(eye),
        sigma2_eye_(sigma_eye * sigma_eye),
        pixel_filter_(std::move(pixel_filter)),
        rgbw_(resolution) {}
  ~RadianceAccumulator() override = default;

  // Resolves the accumulated radiance to the |image| at its given resolution.
  void Resolve(image::Image3f* rgb) const;

  // SampleAccumulator implementation.
  void Add(const RayBundle& bundle,
           absl::Span<const RayBundle::RayIntersectionIndex> solid_samples,
           absl::Span<const int> freespace_rays) override;

 private:
  // Specifies the geometry on which we are representing the surface-lightfield.
  const Frame* frame_;

  // The eye point specifying the direction to resample radiance.
  const ion::math::Point3f eye_;

  // Specifies the falloff of the gaussian kernel used to resample in the
  // direction of 'eye_'.
  const float sigma2_eye_;

  // The filter used to reconstruct pixel color values from samples.
  std::shared_ptr<image::Filter> pixel_filter_;

  // Stores RGB + Weight for additive-blending of all samples.
  image::Image4d rgbw_;
};

}  // namespace baker
}  // namespace seurat

#endif  // VR_SEURAT_BAKER_FRAMEWORK_RADIANCE_ACCUMULATOR_H_
