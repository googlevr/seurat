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

#include "seurat/api/ldi_exporter.h"

#include <cmath>
#include <numeric>

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/image/deep_exr_image_reader.h"

namespace seurat {
namespace api {
namespace {

using ion::math::Vector2i;

TEST(LdiExrIo, WriteDeepExr) {
  Vector2i size(3, 2);
  std::vector<int> sample_counts;
  std::vector<float> colors;
  std::vector<float> depths;

  for (int y = 0; y < size[1]; ++y) {
    for (int x = 0; x < size[0]; ++x) {
      const int sample_count = y * size[0] + x;
      sample_counts.push_back(sample_count);
      for (int i = 0; i < sample_count; ++i) {
        const float value = i / static_cast<float>(sample_count);
        for (int c = 0; c < 4; ++c) {
          colors.push_back(value);
        }
        depths.push_back(value);
      }
    }
  }

  // Write the LDI contents to a temp string
  std::string serialized_exr;
  base::Status status = WriteDeepExr(size[0], size[1], sample_counts, colors,
                                     depths, &serialized_exr);
  EXPECT_TRUE(status.ok());

  // Read the LDI back from the string.
  std::unique_ptr<image::DeepExrImageReader> reader;
  EXPECT_TRUE(image::DeepExrImageReader::Create(serialized_exr, &reader).ok());
  EXPECT_EQ(size, Vector2i(reader->GetImageWidth(), reader->GetImageHeight()));
  EXPECT_EQ(std::accumulate(sample_counts.begin(), sample_counts.end(), 0),
            reader->GetSampleCount());
  EXPECT_EQ(sample_counts, reader->GetSampleCounts());

  std::vector<float> colors_read(reader->GetSampleCount() * 4);
  std::vector<float> depths_read(reader->GetSampleCount());
  reader->AddChannel({"R", &colors_read[0], 4});
  reader->AddChannel({"G", &colors_read[1], 4});
  reader->AddChannel({"B", &colors_read[2], 4});
  reader->AddChannel({"A", &colors_read[3], 4});
  reader->AddChannel({"Z", &depths_read[0], 1});
  EXPECT_TRUE(reader->Read().ok());

  const float kEpsilon = 1e-3;
  EXPECT_EQ(colors.size(), colors_read.size());
  for (size_t i = 0; i < std::min(colors.size(), colors_read.size()); ++i) {
    EXPECT_NEAR(colors[i], colors_read[i], kEpsilon) << i;
  }
  EXPECT_EQ(depths.size(), depths_read.size());
  for (size_t i = 0; i < std::min(depths.size(), depths_read.size()); ++i) {
    EXPECT_NEAR(depths[i], depths_read[i], kEpsilon) << i;
  }
}

}  // namespace
}  // namespace api
}  // namespace seurat
