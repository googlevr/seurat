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

#include "seurat/image/color_processor.h"

#include <array>
#include <memory>

#include "gtest/gtest.h"
#include "seurat/base/camera.h"
#include "seurat/image/image.h"

namespace seurat {
namespace image {
namespace {

using base::Color4f;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector3f;

TEST(ColorProcessors, ToneMapperGammaOneLeavesValuesUnchanged) {
  const float kGamma = 1.0f;
  const GammaToneMapper tone_mapper(kGamma);
  Image4f in_image(2, 2, Color4f::Zero());
  // None of these values require clamping, to simplify the identity test.
  in_image.At(0, 0) = Color4f(0.0f, 0.0f, 0.87f, 1.0f);
  in_image.At(1, 0) = Color4f(0.0f, 1.0f, 0.87f, 0.5f);
  in_image.At(0, 1) = Color4f(1.0f, 1.0f, 0.0f, 0.0f);
  in_image.At(1, 1) = Color4f(0.25f, 0.1f, 1.0f, 1.0f);

  Image4f image = in_image;
  tone_mapper.ProcessColors(
      absl::Span<Color4f>(image.Data(), image.Width() * image.Height()));

  for (int y = 0; y < image.Height(); ++y) {
    for (int x = 0; x < image.Width(); ++x) {
      EXPECT_EQ(in_image.At(x, y), image.At(x, y)) << Point2i(x, y);
    }
  }
}

TEST(ColorProcessors, ToneMapperGammaTwoTakesSqrtOfValues) {
  const float kGamma = 2.0f;
  const GammaToneMapper tone_mapper(kGamma);
  Image4f in_image(2, 2, Color4f::Zero());
  // None of these values require clamping, to simplify the square-root test.
  in_image.At(0, 0) = Color4f(0.01f, 0.99f, 0.87f, 1.0f);
  in_image.At(1, 0) = Color4f(0.0f, 1.0f, 0.87f, 0.5f);
  in_image.At(0, 1) = Color4f(1.0f, 1.0f, 0.5f, 0.0f);
  in_image.At(1, 1) = Color4f(0.25f, 0.1f, 1.0f, 1.0f);

  Image4f image = in_image;
  tone_mapper.ProcessColors(
      absl::Span<Color4f>(image.Data(), image.Width() * image.Height()));

  for (int y = 0; y < image.Height(); ++y) {
    for (int x = 0; x < image.Width(); ++x) {
      const Color4f in_pixel = in_image.At(x, y);
      Color4f color_sqrt(std::sqrt(in_pixel[0]), std::sqrt(in_pixel[1]),
                         std::sqrt(in_pixel[2]), in_pixel[3]);
      EXPECT_EQ(image.At(x, y), color_sqrt) << Point2i(x, y);
    }
  }
}

TEST(ColorProcessors, PremultipliedAlphaConverter) {
  PremultipliedAlphaConverter pma;
  std::array<Color4f, 1> test_colors = {{Color4f(1.0f, 0.0f, 2.0f, 0.5f)}};
  pma.ProcessColors(absl::MakeSpan(test_colors));
  EXPECT_EQ(Color4f(0.5f, 0.0f, 1.0f, 0.5f), test_colors[0]);
}

TEST(ColorProcessors, ColorProcessorPipeline) {
  // For this test, convert to premultiplied alpha first, and then tone-map.
  std::vector<std::unique_ptr<ColorProcessor>> sequence;
  sequence.emplace_back(new PremultipliedAlphaConverter);
  sequence.emplace_back(new GammaToneMapper(2.2f));
  ColorProcessorPipeline pipeline(std::move(sequence));

  // Premultiplying by alpha & then tone-mapping should result in identical RGB
  // values for both pixels.
  std::array<Color4f, 2> test_colors = {
      {Color4f(1.0f, 0.0f, 0.4f, 0.5f), Color4f(0.5f, 0.0f, 0.2f, 1.0f)}};
  pipeline.ProcessColors(absl::MakeSpan(test_colors));
  EXPECT_EQ(test_colors[0][0], test_colors[1][0]);
  EXPECT_EQ(test_colors[0][1], test_colors[1][1]);
  EXPECT_EQ(test_colors[0][2], test_colors[1][2]);
}

}  // namespace
}  // namespace image
}  // namespace seurat
