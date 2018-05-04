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

#ifndef VR_SEURAT_IMAGE_COLOR_PROCESSOR_H_
#define VR_SEURAT_IMAGE_COLOR_PROCESSOR_H_

#include "absl/types/span.h"
#include "seurat/base/color.h"

namespace seurat {
namespace image {

// ColorProcessor specifies an interface for transforming sets of colors.
//
// Implementations may perform operations such as tone-mapping and/or
// transformation to premultiplied alpha.
class ColorProcessor {
 public:
  virtual ~ColorProcessor() = default;

  // Processes the color values in |color_set| in place.
  virtual void ProcessColors(absl::Span<base::Color4f> color_set) const = 0;

 protected:
  ColorProcessor() = default;
};

// GammaToneMapper tonemaps sets of colors via a gamma exponent curve. The
// process expects display-referred (i.e. sRGB linear or gamma space),
// possibly-HDR colors as input, and produces the same kind of data, with a
// power distribution depending on the configured gamma-correction value.
class GammaToneMapper : public ColorProcessor {
 public:
  // Constructs a tone mapper with the given exponent for |gamma| correction.
  explicit GammaToneMapper(float gamma) : gamma_(gamma) {}

  ~GammaToneMapper() override = default;

  void ProcessColors(absl::Span<base::Color4f> color_set) const override;

 private:
  // Exponent for gamma correction.
  const float gamma_;
};

// Transforms independent rgba into premultiplied-alpha.
class PremultipliedAlphaConverter : public ColorProcessor {
 public:
  PremultipliedAlphaConverter() = default;

  ~PremultipliedAlphaConverter() override = default;

  void ProcessColors(absl::Span<base::Color4f> color_set) const override;
};

// Applies a sequence of ColorProcessors in order.
class ColorProcessorPipeline : public ColorProcessor {
 public:
  explicit ColorProcessorPipeline(
      std::vector<std::unique_ptr<ColorProcessor>> sequence)
      : sequence_(std::move(sequence)) {}

  ~ColorProcessorPipeline() override = default;

  void ProcessColors(absl::Span<base::Color4f> color_set) const override;

 private:
  const std::vector<std::unique_ptr<ColorProcessor>> sequence_;
};

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_COLOR_PROCESSOR_H_
