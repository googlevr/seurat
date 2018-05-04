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

#include "seurat/image/ion_image_reader.h"
#include "seurat/image/ion_image_writer.h"

#include <cmath>   // NOLINT
#include <memory>  // NOLINT
#include <string>  // NOLINT

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
      const float value = ((x + y * size[0]) % 252) / 255.0f;
      image.At(x, y) = Color4f(value,                  //
                               1.0f / 255.0f + value,  //
                               2.0f / 255.0f + value,  //
                               3.0f / 255.0f + value);
    }
  }
  return image;
}

base::Status ReadPng(const std::string& serialized_png, Image4f* image) {
  // Read an image back from the string.
  std::unique_ptr<IonImageReader> reader;
  SEURAT_RETURN_IF_ERROR(IonImageReader::Create(serialized_png, &reader));
  image->Resize(reader->GetImageWidth(), reader->GetImageHeight());
  const int dimension = Image4f::ElementType::kDimension;
  reader->AddChannel({"R", &image->data()[0][0], dimension});
  reader->AddChannel({"G", &image->data()[0][1], dimension});
  reader->AddChannel({"B", &image->data()[0][2], dimension});
  reader->AddChannel({"A", &image->data()[0][3], dimension});
  SEURAT_RETURN_IF_ERROR(reader->Read());
  return base::OkStatus();
}

base::Status WritePng(const Image4f& image, std::string* serialized_png) {
  return WriteImagePngToString(image, serialized_png);
}

TEST(IonIo, WriteEmptyPngFails) {
  Vector2i size(0, 0);
  Image4f image(size);

  // Try writing the image to a temp string.
  std::string serialized_png;
  base::Status status = WritePng(image, &serialized_png);
  EXPECT_FALSE(status.ok());
  EXPECT_EQ("Image size (0, 0) is not positive", status.error_message());
}

TEST(IonIo, WriteAndReadPng) {
  // Create a test image and write it to a temp string.
  const Image4f image = CreateTestImage({3, 2});
  std::string serialized_png;
  EXPECT_TRUE(WritePng(image, &serialized_png).ok());

  // Read the image back from the string.
  Image4f image_copy;
  EXPECT_TRUE(ReadPng(serialized_png, &image_copy).ok());

  EXPECT_EQ(image.GetSize(), image_copy.GetSize());
  for (int y = 0; y < image.Height(); ++y) {
    for (int x = 0; x < image.Width(); ++x) {
      EXPECT_TRUE(
          ColorsAlmostEqual(image.At(x, y), image_copy.At(x, y), 1e-6f));
    }
  }
}

TEST(IonIo, ReadSpecialChannels) {
  // Create a test image and write it to a temp string.
  const Image4f image = CreateTestImage({3, 2});
  std::string serialized_png;
  EXPECT_TRUE(WritePng(image, &serialized_png).ok());

  // Read the image back from the string, but swizzle components and also use
  // some of the reserved channels.
  std::unique_ptr<IonImageReader> reader;
  EXPECT_TRUE(IonImageReader::Create(serialized_png, &reader).ok());
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

  EXPECT_EQ(modified_image.GetSize(), image_copy.GetSize());
  for (int y = 0; y < modified_image.Height(); ++y) {
    for (int x = 0; x < modified_image.Width(); ++x) {
      EXPECT_TRUE(ColorsAlmostEqual(modified_image.At(x, y),
                                    image_copy.At(x, y), 1e-6f));
    }
  }
}

TEST(IonIo, NonExistentChannels) {
  // Create a test image and write it to a temp string.
  const Image4f image = CreateTestImage({3, 2});
  std::string serialized_png;
  EXPECT_TRUE(WritePng(image, &serialized_png).ok());

  // Read the image back from the string, but specify a channel that does not
  // exist.
  std::unique_ptr<IonImageReader> reader;
  EXPECT_TRUE(IonImageReader::Create(serialized_png, &reader).ok());
  Image4f image_copy(reader->GetImageWidth(), reader->GetImageHeight());
  const int dimension = Image4f::ElementType::kDimension;
  reader->AddChannel({"R", &image_copy.data()[0][0], dimension});
  reader->AddChannel({"foo", &image_copy.data()[0][1], dimension});
  reader->AddChannel({"B", &image_copy.data()[0][2], dimension});
  reader->AddChannel({"A", &image_copy.data()[0][3], dimension});
  base::Status status = reader->Read();
  EXPECT_FALSE(status.ok());
  EXPECT_EQ("Invalid channel name: foo", status.error_message());
}

TEST(ImageLoaderTest, LoadPngColorImage) {
  std::shared_ptr<base::FileSystem> file_system =
      std::make_shared<base::FileSystem>(testing::GetTestSrcdir());
  const char kPngFilename[] = "com_google_seurat/seurat/image/testdata/color.png";
  std::string serialized_png;
  EXPECT_TRUE(file_system->GetContents(kPngFilename, &serialized_png).ok());

  std::unique_ptr<IonImageReader> reader;
  EXPECT_TRUE(IonImageReader::Create(serialized_png, &reader).ok());
  Image3f image(reader->GetImageWidth(), reader->GetImageHeight());
  const int dimension = Image3f::ElementType::kDimension;
  reader->AddChannel({"R", &image.data()[0][0], dimension});
  reader->AddChannel({"G", &image.data()[0][1], dimension});
  reader->AddChannel({"B", &image.data()[0][2], dimension});
  EXPECT_TRUE(reader->Read().ok());

  // Pixel at (0, 511) has a greyish color of (104, 97, 91).
  const Color3f kExpectedColor_0_511 = Color3ui8(104, 97, 91).AsColorF();
  EXPECT_EQ(kExpectedColor_0_511, image.At(0, 511));
  // Pixel at (256, 256) to be red.
  const Color3f kExpectedColor_256_256 = Color3f(1.0f, 0.0f, 0.0f);
  EXPECT_EQ(kExpectedColor_256_256, image.At(256, 256));
}

}  // namespace
}  // namespace image
}  // namespace seurat
