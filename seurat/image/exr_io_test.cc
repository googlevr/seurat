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

#include "seurat/image/exr_image_reader.h"
#include "seurat/image/exr_image_writer.h"

#include <cmath>
#include <memory>
#include <string>

#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "gtest/gtest.h"
#include "seurat/base/array2d_util.h"
#include "seurat/base/file_system.h"
#include "seurat/base/status.h"
#include "seurat/image/image.h"
#include "seurat/testing/ion_test_utils.h"
#include "seurat/testing/test_flags.h"

namespace seurat {
namespace image {
namespace {

using base::Color1f;
using base::Color3f;
using base::Color3ui8;
using base::Color4f;
using ion::math::Point2i;
using ion::math::Vector2i;

Image4f CreateTestImage(const Vector2i& size) {
  Image4f image(size[0], size[1]);
  for (int y = 0; y < size[1]; ++y) {
    for (int x = 0; x < size[0]; ++x) {
      const float value = 1.0f / static_cast<float>(x + y * size[0]);
      image.At(x, y) = Color4f(value, value + 1.0f, value + 2.0f, value + 3.0f);
    }
  }
  return image;
}

// Comparator for color values, taking into account float->half precision loss.
template <typename ColorType>
bool HalfColorAlmostEqual(const ColorType& lhs, const ColorType& rhs) {
  const float kEpsilon = 8e-4f;
  for (int i = 0; i < ColorType::kDimension; ++i) {
    if (std::abs(lhs[i] - rhs[i]) > kEpsilon) {
      return false;
    }
  }
  return true;
}

base::Status ReadExr(const std::string& serialized_exr, Image4f* image) {
  // Read an image back from the string.
  std::unique_ptr<ExrImageReader> reader;
  SEURAT_RETURN_IF_ERROR(ExrImageReader::Create(serialized_exr, &reader));
  image->Resize(reader->GetImageWidth(), reader->GetImageHeight());
  const int dimension = Image4f::ElementType::kDimension;
  reader->AddChannel({"R", &image->data()[0][0], dimension});
  reader->AddChannel({"G", &image->data()[0][1], dimension});
  reader->AddChannel({"B", &image->data()[0][2], dimension});
  reader->AddChannel({"A", &image->data()[0][3], dimension});
  SEURAT_RETURN_IF_ERROR(reader->Read());
  return base::OkStatus();
}

base::Status WriteExr(const Image4f& image,
                      ExrImageWriter::ValueType color_type,
                      std::string* serialized_exr) {
  // Write the Image to a string.
  std::unique_ptr<ExrImageWriter> writer;
  SEURAT_RETURN_IF_ERROR(
      ExrImageWriter::Create(image.Width(), image.Height(), &writer));
  const int dimension = Image4f::ElementType::kDimension;
  writer->AddChannel({"R", &image.data()[0][0], dimension, color_type});
  writer->AddChannel({"G", &image.data()[0][1], dimension, color_type});
  writer->AddChannel({"B", &image.data()[0][2], dimension, color_type});
  writer->AddChannel({"A", &image.data()[0][3], dimension, color_type});
  SEURAT_RETURN_IF_ERROR(writer->Write(serialized_exr));
  return base::OkStatus();
}

TEST(ExrIo, WriteEmptyExrFails) {
  Vector2i size(0, 0);
  Image4f image(size);

  // Try writing the image to a temp string.
  std::string serialized_exr;
  base::Status status =
      WriteExr(image, ExrImageWriter::ValueType::kFloat, &serialized_exr);
  EXPECT_FALSE(status.ok());
  EXPECT_EQ("EXR size (0, 0) is not positive", status.error_message());
}

TEST(ExrIo, ReadInvalidExrFails) {
  std::string empty_string;
  Image4f image;
  base::Status status = ReadExr(empty_string, &image);
  EXPECT_FALSE(status.ok());
}

TEST(ExrIo, WriteAndReadExr) {
  // Create a test image and write it to a temp string.
  const Image4f image = CreateTestImage({3, 2});
  std::string serialized_exr;
  EXPECT_TRUE(
      WriteExr(image, ExrImageWriter::ValueType::kFloat, &serialized_exr).ok());

  // Read the image back from the string.
  Image4f image_copy;
  EXPECT_TRUE(ReadExr(serialized_exr, &image_copy).ok());

  EXPECT_EQ(image, image_copy);
}

TEST(ExrIo, ReadSpecialChannels) {
  // Create a test image and write it to a temp string.
  const Image4f image = CreateTestImage({3, 2});
  std::string serialized_exr;
  EXPECT_TRUE(
      WriteExr(image, ExrImageWriter::ValueType::kFloat, &serialized_exr).ok());

  // Read the image back from the string, but swizzle components and also use
  // some of the reserved channels.
  std::unique_ptr<ExrImageReader> reader;
  EXPECT_TRUE(ExrImageReader::Create(serialized_exr, &reader).ok());
  Image4f image_copy(reader->GetImageWidth(), reader->GetImageHeight());
  const int dimension = Image4f::ElementType::kDimension;
  reader->AddChannel({"G", &image_copy.data()[0][0], dimension});
  reader->AddChannel({"B", &image_copy.data()[0][1], dimension});
  reader->AddChannel({"CONSTANT_ZERO", &image_copy.data()[0][2], dimension});
  reader->AddChannel({"CONSTANT_ONE", &image_copy.data()[0][3], dimension});
  EXPECT_TRUE(reader->Read().ok());

  Image4f modified_image(image.Width(), image.Height());
  base::TransformArray(image, &modified_image, [](const Color4f& color) {
    return Color4f(color[1], color[2], 0.0f, 1.0f);
  });
  EXPECT_EQ(modified_image, image_copy);
}

TEST(ExrIo, NonExistentChannels) {
  // Create a test image and write it to a temp string.
  const Image4f image = CreateTestImage({3, 2});
  std::string serialized_exr;
  EXPECT_TRUE(
      WriteExr(image, ExrImageWriter::ValueType::kFloat, &serialized_exr).ok());

  // Read the image back from the string, but specify a channel that does not
  // exist.
  std::unique_ptr<ExrImageReader> reader;
  EXPECT_TRUE(ExrImageReader::Create(serialized_exr, &reader).ok());
  Image4f image_copy(reader->GetImageWidth(), reader->GetImageHeight());
  const int dimension = Image4f::ElementType::kDimension;
  reader->AddChannel({"R", &image_copy.data()[0][0], dimension});
  reader->AddChannel({"foo", &image_copy.data()[0][1], dimension});
  reader->AddChannel({"B", &image_copy.data()[0][2], dimension});
  reader->AddChannel({"A", &image_copy.data()[0][3], dimension});
  base::Status status = reader->Read();
  EXPECT_FALSE(status.ok());
  EXPECT_EQ("EXR channel \"foo\" missing. Existing channels: A, B, G, R",
            status.error_message());
}

TEST(ExrIo, HalfPrecisionConversion) {
  // Create a test image and write it to a temp string.
  const Image4f image = CreateTestImage({3, 2});
  std::string serialized_exr;
  EXPECT_TRUE(
      WriteExr(image, ExrImageWriter::ValueType::kHalf, &serialized_exr).ok());

  // Read the image back from the string.
  Image4f image_copy;
  EXPECT_TRUE(ReadExr(serialized_exr, &image_copy).ok());

  EXPECT_EQ(image.GetSize(), image_copy.GetSize());
  for (int y = 0; y < image.Height(); ++y) {
    for (int x = 0; x < image.Width(); ++x) {
      EXPECT_TRUE(HalfColorAlmostEqual(image.At(x, y), image_copy.At(x, y)));
    }
  }
}

TEST(ImageLoaderTest, LoadExrColorImage) {
  std::shared_ptr<base::FileSystem> file_system =
      std::make_shared<base::FileSystem>(testing::GetTestSrcdir());
  const char kExrFilename[] = "com_google_seurat/seurat/image/testdata/color.exr";
  std::string serialized_exr;
  EXPECT_TRUE(file_system->GetContents(kExrFilename, &serialized_exr).ok());

  std::unique_ptr<ExrImageReader> reader;
  EXPECT_TRUE(ExrImageReader::Create(serialized_exr, &reader).ok());
  Image3f image(reader->GetImageWidth(), reader->GetImageHeight());
  const int dimension = Image3f::ElementType::kDimension;
  reader->AddChannel({"R", &image.data()[0][0], dimension});
  reader->AddChannel({"G", &image.data()[0][1], dimension});
  reader->AddChannel({"B", &image.data()[0][2], dimension});
  EXPECT_TRUE(reader->Read().ok());
  base::FlipArrayVertical(&image);

  // Pixel at (0, 0) has a greyish color.
  const Color3f kExpectedColor_0_0(0.410156f, 0.381836f, 0.357910f);
  constexpr float kColorTolerance = 1.0f / 255;
  EXPECT_VECTOR_NEAR(kExpectedColor_0_0, image.At(0, 0), kColorTolerance);
  // Pixel at (245, 217) is off-red.
  const Color3f kExpectedColor_245_217 =
      Color3f(0.978516f, 0.178101f, 0.184814f);
  EXPECT_VECTOR_NEAR(kExpectedColor_245_217, image.At(245, 217),
                     kColorTolerance);
}

TEST(ImageLoaderTest, LoadExrDepthImage) {
  std::shared_ptr<base::FileSystem> file_system =
      std::make_shared<base::FileSystem>(testing::GetTestSrcdir());
  const char kExrFilename[] = "com_google_seurat/seurat/image/testdata/depth.exr";
  std::string serialized_exr;
  EXPECT_TRUE(file_system->GetContents(kExrFilename, &serialized_exr).ok());

  std::unique_ptr<ExrImageReader> reader;
  EXPECT_TRUE(ExrImageReader::Create(serialized_exr, &reader).ok());
  Image1f image(reader->GetImageWidth(), reader->GetImageHeight());
  const int dimension = Image1f::ElementType::kDimension;
  reader->AddChannel({"R", &image.data()[0][0], dimension});
  EXPECT_TRUE(reader->Read().ok());
  base::FlipArrayVertical(&image);

  // Depth at (275, 275) is ~15.0 units away.
  const Color1f kExpectedColor_274_274 = Color1f(15.0f);
  EXPECT_VECTOR_FLOAT_EQ(kExpectedColor_274_274, image.At(274, 274));

  // Depth at (245, 217) is ~10.0 unit away.
  const Color1f kExpectedColor_245_217 = Color1f(10.0f);
  EXPECT_VECTOR_FLOAT_EQ(kExpectedColor_245_217, image.At(245, 217));

  // Depth at (303, 232) is ~11.0 units away.
  const Color1f kExpectedColor_303_232 = Color1f(11.0f);
  EXPECT_VECTOR_FLOAT_EQ(kExpectedColor_303_232, image.At(303, 232));
}

}  // namespace
}  // namespace image
}  // namespace seurat
