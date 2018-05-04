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

#include "seurat/image/deep_exr_image_reader.h"
#include "seurat/image/deep_exr_image_writer.h"

#include <cmath>
#include <memory>
#include <string>

#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "gtest/gtest.h"
#include "seurat/base/status.h"
#include "seurat/image/ldi.h"
#include "seurat/image/ldi_test_utils.h"
#include "seurat/testing/ion_test_utils.h"
#include "seurat/testing/test_flags.h"

namespace seurat {
namespace image {
namespace {

using base::Color4f;
using ion::math::Point2i;
using ion::math::Vector2i;

Ldi4f CreateTestLdi(const Vector2i& size) {
  std::vector<int> sample_counts;
  std::vector<Color4f> colors;
  std::vector<float> depths;
  for (int y = 0; y < size[1]; ++y) {
    for (int x = 0; x < size[0]; ++x) {
      const int sample_count = y * size[0] + x;
      sample_counts.push_back(sample_count);
      for (int i = 0; i < sample_count; ++i) {
        const float value = i / static_cast<float>(sample_count);
        colors.push_back(Color4f(value, value, value, value));
        depths.push_back(value);
      }
    }
  }
  return Ldi4f(size, sample_counts, colors, depths);
}

// Comparator for color values, taking into account float->half precision loss.
template <typename ColorType>
bool HalfColorAlmostEqual(const ColorType& lhs, const ColorType& rhs) {
  const float kEpsilon = 2e-4f;
  for (int i = 0; i < ColorType::kDimension; ++i) {
    if (std::abs(lhs[i] - rhs[i]) > kEpsilon) {
      return false;
    }
  }
  return true;
}

// Comparator for depth values, taking into account float->half precision loss.
template <typename T>
bool HalfDepthAlmostEqual(const T& lhs, const T& rhs) {
  const float kEpsilon = 2e-4f;
  return (std::abs(lhs - rhs) <= kEpsilon);
}

base::Status ReadExr(const std::string& serialized_exr, Ldi4f* ldi) {
  // Read an LDI back from the string.
  std::unique_ptr<DeepExrImageReader> reader;
  SEURAT_RETURN_IF_ERROR(
      DeepExrImageReader::Create(serialized_exr, &reader));
  std::vector<Color4f> colors(reader->GetSampleCount());
  std::vector<float> depths(reader->GetSampleCount());
  const auto sample_counts = reader->GetSampleCounts();
  reader->AddChannel({"R", &colors.data()[0][0], 4});
  reader->AddChannel({"G", &colors.data()[0][1], 4});
  reader->AddChannel({"B", &colors.data()[0][2], 4});
  reader->AddChannel({"A", &colors.data()[0][3], 4});
  reader->AddChannel({"Z", &depths.data()[0], 1});
  SEURAT_RETURN_IF_ERROR(reader->Read());

  *ldi = Ldi4f({reader->GetImageWidth(), reader->GetImageHeight()},
               {sample_counts.begin(), sample_counts.end()}, std::move(colors),
               std::move(depths));
  return base::OkStatus();
}

base::Status WriteExr(const Ldi4f& ldi,
                      DeepExrImageWriter::ValueType color_type,
                      DeepExrImageWriter::ValueType depth_type,
                      std::string* serialized_exr) {
  // Write the LDI to a string.
  std::vector<int> sample_counts;
  sample_counts.reserve(ldi.GetWidth() * ldi.GetHeight());
  for (int y = 0; y < ldi.GetHeight(); ++y) {
    for (int x = 0; x < ldi.GetWidth(); ++x) {
      sample_counts.push_back(ldi.GetSampleCount({x, y}));
    }
  }
  std::unique_ptr<DeepExrImageWriter> writer;
  SEURAT_RETURN_IF_ERROR(DeepExrImageWriter::Create(
      ldi.GetSize()[0], ldi.GetSize()[1], sample_counts, &writer));
  writer->AddChannel({"R", &ldi.GetColors()[0][0], 4, color_type});
  writer->AddChannel({"G", &ldi.GetColors()[0][1], 4, color_type});
  writer->AddChannel({"B", &ldi.GetColors()[0][2], 4, color_type});
  writer->AddChannel({"A", &ldi.GetColors()[0][3], 4, color_type});
  writer->AddChannel({"Z", &ldi.GetDepths()[0], 1, depth_type});
  SEURAT_RETURN_IF_ERROR(writer->Write(serialized_exr));
  return base::OkStatus();
}

TEST(DeepExrIo, WriteAndReadEmptyExr) {
  Vector2i size(3, 2);
  std::vector<int> sample_counts{{0, 0, 0, 0, 0, 0}};
  std::vector<Color4f> colors;
  std::vector<float> depths;
  Ldi4f ldi(size, sample_counts, colors, depths);

  // Write the LDI to a temp string
  std::string serialized_exr;
  base::Status status =
      WriteExr(ldi, DeepExrImageWriter::ValueType::kFloat,
               DeepExrImageWriter::ValueType::kFloat, &serialized_exr);
  EXPECT_TRUE(status.ok());

  // Read the LDI back from the string.
  Ldi4f ldi_copy;
  EXPECT_TRUE(ReadExr(serialized_exr, &ldi_copy).ok());

  EXPECT_TRUE(LdiEquals(ldi, ldi_copy));
}

TEST(DeepExrIo, WriteAndReadExrWithEmptyPixelAtEndOfScanline) {
  // This particular configuration was triggering a bug in the past, when the
  // sample count slice was stubbed out by a dummy in DeepExrImageReader::Read,
  // incorrectly assuming that OpenEXR would use the data from the previous call
  // in DeepExrImageReader::Init.
  Vector2i size(3, 2);
  std::vector<int> sample_counts{{1, 1, 1, 1, 1, 0}};
  std::vector<Color4f> colors{{
      Color4f(1.0f, 0.0f, 0.0f, 1.0f),
      Color4f(1.0f, 0.0f, 0.0f, 1.0f),
      Color4f(1.0f, 0.0f, 0.0f, 1.0f),
      Color4f(1.0f, 0.0f, 0.0f, 1.0f),
      Color4f(1.0f, 0.0f, 0.0f, 1.0f),
  }};
  std::vector<float> depths{{0.5f, 0.5f, 0.5f, 0.5f, 0.5f}};
  Ldi4f ldi(size, sample_counts, colors, depths);

  // Write the LDI to a temp string
  std::string serialized_exr;
  base::Status status =
      WriteExr(ldi, DeepExrImageWriter::ValueType::kFloat,
               DeepExrImageWriter::ValueType::kFloat, &serialized_exr);
  EXPECT_TRUE(status.ok());

  // Read the LDI back from the string.
  Ldi4f ldi_copy;
  EXPECT_TRUE(ReadExr(serialized_exr, &ldi_copy).ok());

  EXPECT_TRUE(LdiEquals(ldi, ldi_copy));
}

TEST(DeepExrIo, WriteAndReadExr) {
  // Create a test LDI and write it to a temp string.
  const Ldi4f ldi = CreateTestLdi({3, 2});
  std::string serialized_exr;
  base::Status status =
      WriteExr(ldi, DeepExrImageWriter::ValueType::kFloat,
               DeepExrImageWriter::ValueType::kFloat, &serialized_exr);
  EXPECT_TRUE(status.ok());

  // Read the LDI back from the string.
  Ldi4f ldi_copy;
  EXPECT_TRUE(ReadExr(serialized_exr, &ldi_copy).ok());

  EXPECT_TRUE(LdiEquals(ldi, ldi_copy));
}

TEST(DeepExrIo, ReadSpecialChannels) {
  // Create a test LDI and write it to a temp string.
  const Ldi4f ldi = CreateTestLdi({3, 2});
  std::string serialized_exr;
  base::Status status =
      WriteExr(ldi, DeepExrImageWriter::ValueType::kFloat,
               DeepExrImageWriter::ValueType::kFloat, &serialized_exr);
  EXPECT_TRUE(status.ok());

  // Read the LDI back from the string, but swizzle components and also use some
  // of the reserved channels.
  std::unique_ptr<DeepExrImageReader> reader;
  EXPECT_TRUE(DeepExrImageReader::Create(serialized_exr, &reader).ok());
  std::vector<Color4f> read_colors(reader->GetSampleCount());
  std::vector<float> read_depths(reader->GetSampleCount());
  const auto sample_counts = reader->GetSampleCounts();
  reader->AddChannel({"G", &read_colors[0][0], 4});
  reader->AddChannel({"B", &read_colors[0][1], 4});
  reader->AddChannel({"CONSTANT_ZERO", &read_colors[0][2], 4});
  reader->AddChannel({"CONSTANT_ONE", &read_colors[0][3], 4});
  reader->AddChannel({"Z", &read_depths[0], 1});
  EXPECT_TRUE(reader->Read().ok());
  Ldi4f ldi_copy = Ldi4f({reader->GetImageWidth(), reader->GetImageHeight()},
                         {sample_counts.begin(), sample_counts.end()},
                         std::move(read_colors), std::move(read_depths));

  Ldi4f modified_ldi = ldi;
  for (auto& color : modified_ldi.GetMutableColors()) {
    color = Color4f(color[1], color[2], 0.0f, 1.0f);
  }
  EXPECT_TRUE(LdiEquals(modified_ldi, ldi_copy, HalfColorAlmostEqual<Color4f>));
}

TEST(DeepExrIo, NonExistentChannels) {
  // Create a test LDI and write it to a temp string.
  const Ldi4f ldi = CreateTestLdi({3, 2});
  std::string serialized_exr;
  base::Status write_status =
      WriteExr(ldi, DeepExrImageWriter::ValueType::kFloat,
               DeepExrImageWriter::ValueType::kFloat, &serialized_exr);
  EXPECT_TRUE(write_status.ok());

  // Read the LDI back from the string, but specify a channel that does not
  // exist.
  std::unique_ptr<DeepExrImageReader> reader;
  EXPECT_TRUE(DeepExrImageReader::Create(serialized_exr, &reader).ok());
  std::vector<Color4f> read_colors(reader->GetSampleCount());
  std::vector<float> read_depths(reader->GetSampleCount());
  reader->AddChannel({"R", &read_colors[0][0], 4});
  reader->AddChannel({"foo", &read_colors[0][1], 4});
  reader->AddChannel({"B", &read_colors[0][2], 4});
  reader->AddChannel({"A", &read_colors[0][3], 4});
  reader->AddChannel({"Z", &read_depths[0], 1});
  base::Status read_status = reader->Read();
  EXPECT_FALSE(read_status.ok());
  EXPECT_EQ("EXR channel \"foo\" missing", read_status.error_message());
}

TEST(DeepExrIo, HalfPrecisionConversion) {
  // Create a test LDI and write it to a temp string.
  const Ldi4f ldi = CreateTestLdi({3, 2});
  std::string serialized_exr;
  base::Status status =
      WriteExr(ldi, DeepExrImageWriter::ValueType::kHalf,
               DeepExrImageWriter::ValueType::kHalf, &serialized_exr);
  EXPECT_TRUE(status.ok());

  // Read the LDI back from the string.
  Ldi4f ldi_copy;
  EXPECT_TRUE(ReadExr(serialized_exr, &ldi_copy).ok());

  EXPECT_TRUE(LdiEquals(ldi, ldi_copy, HalfColorAlmostEqual<Color4f>,
                        HalfDepthAlmostEqual<float>));
}

}  // namespace
}  // namespace image
}  // namespace seurat
