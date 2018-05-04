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

#include "seurat/image/image_util.h"

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/color.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace image {
namespace {

using base::Color1f;
using base::Color3f;
using base::Color4f;
using ion::math::Point2i;
using ion::math::Range2i;
using ion::math::Vector2i;

constexpr float kEpsilon = 1.0e-5f;

TEST(ImageUtilTest, ConvertSeuratImageToIonImage_Image4f) {
  Image4f simage(16, 8);
  for (int y = 0; y < 8; ++y) {
    for (int x = 0; x < 16; ++x) {
      simage.At(x, y) = {0.0f, 0.5f, 1.0f, 0.5f};
    }
  }
  const ion::gfx::ImagePtr ion_image = ConvertSeuratImageToIonImage(simage);
  EXPECT_EQ(16, ion_image->GetWidth());
  EXPECT_EQ(8, ion_image->GetHeight());
  EXPECT_EQ(ion::gfx::Image::Format::kRgba8888, ion_image->GetFormat());
  const uint8* ion_image_data = ion_image->GetData()->GetData<uint8>();
  for (int i = 0; i < 16 * 8; ++i) {
    EXPECT_EQ(0, ion_image_data[i * 4 + 0]);
    EXPECT_EQ(128, ion_image_data[i * 4 + 1]);
    EXPECT_EQ(255, ion_image_data[i * 4 + 2]);
    EXPECT_EQ(128, ion_image_data[i * 4 + 3]);
  }
}

TEST(ImageUtilTest, ConvertSeuratImageToIonImage_Image3f) {
  Image3f simage(16, 8);
  for (int y = 0; y < 8; ++y) {
    for (int x = 0; x < 16; ++x) {
      simage.At(x, y) = {0.0f, 0.5f, 1.0f};
    }
  }
  const ion::gfx::ImagePtr ion_image = ConvertSeuratImageToIonImage(simage);
  EXPECT_EQ(16, ion_image->GetWidth());
  EXPECT_EQ(8, ion_image->GetHeight());
  EXPECT_EQ(ion::gfx::Image::Format::kRgb888, ion_image->GetFormat());
  const uint8* ion_image_data = ion_image->GetData()->GetData<uint8>();
  for (int i = 0; i < 16 * 8; ++i) {
    EXPECT_EQ(0, ion_image_data[i * 3 + 0]);
    EXPECT_EQ(128, ion_image_data[i * 3 + 1]);
    EXPECT_EQ(255, ion_image_data[i * 3 + 2]);
  }
}

TEST(ImageUtilTest, ConvertSeuratImageToIonImage_Image1f) {
  Image1f simage(16, 8);
  for (int y = 0; y < 8; ++y) {
    for (int x = 0; x < 16; ++x) {
      simage.At(x, y) = Color1f(0.5f);
    }
  }
  const ion::gfx::ImagePtr ion_image = ConvertSeuratImageToIonImage(simage);
  EXPECT_EQ(16, ion_image->GetWidth());
  EXPECT_EQ(8, ion_image->GetHeight());
  EXPECT_EQ(ion::gfx::Image::Format::kAlpha, ion_image->GetFormat());
  const uint8* ion_image_data = ion_image->GetData()->GetData<uint8>();
  for (int i = 0; i < 16 * 8; ++i) {
    EXPECT_EQ(128, ion_image_data[i]);
  }
}

TEST(ImageUtilTest, ConvertSeuratImageToIonImage_Image4ui8) {
  Image4ui8 simage(16, 8);
  for (int y = 0; y < 8; ++y) {
    for (int x = 0; x < 16; ++x) {
      simage.At(x, y) = Image4ui8::ElementType(x, y, x, y);
    }
  }
  const ion::gfx::ImagePtr ion_image = ConvertSeuratImageToIonImage(simage);
  EXPECT_EQ(16, ion_image->GetWidth());
  EXPECT_EQ(8, ion_image->GetHeight());
  EXPECT_EQ(ion::gfx::Image::Format::kRgba8888, ion_image->GetFormat());
  const uint8* ion_image_data = ion_image->GetData()->GetData<uint8>();
  for (int i = 0; i < 16 * 8; ++i) {
    EXPECT_EQ(i % 16, ion_image_data[i * 4 + 0]);
    EXPECT_EQ(i / 16, ion_image_data[i * 4 + 1]);
    EXPECT_EQ(i % 16, ion_image_data[i * 4 + 2]);
    EXPECT_EQ(i / 16, ion_image_data[i * 4 + 3]);
  }
}

TEST(ImageUtilTest, ConvertIonImageToWImage4b) {
  const int kWidth = 16;
  const int kHeight = 8;
  ion::gfx::ImagePtr ion_image =
      base::CreateImage(ion::gfx::Image::kRgba8888, {kWidth, kHeight});
  const auto data_container = ion_image->GetData();
  base::Color4ui8* data = data_container->GetMutableData<base::Color4ui8>();
  for (int y = 0; y < kHeight; ++y) {
    for (int x = 0; x < kWidth; ++x) {
      data[y * kWidth + x][0] = 0;
      data[y * kWidth + x][1] = 128;
      data[y * kWidth + x][2] = 255;
      data[y * kWidth + x][3] = 128;
    }
  }

  const auto simage = ConvertIonImageToSeuratImage<Image4ui8>(ion_image);
  EXPECT_EQ(kWidth, simage.Width());
  EXPECT_EQ(kHeight, simage.Height());
  for (int y = 0; y < kHeight; ++y) {
    for (int x = 0; x < kWidth; ++x) {
      EXPECT_EQ(0, simage.At(x, y)[0]);
      EXPECT_EQ(128, simage.At(x, y)[1]);
      EXPECT_EQ(255, simage.At(x, y)[2]);
      EXPECT_EQ(128, simage.At(x, y)[3]);
    }
  }
}
}  // namespace
}  // namespace image
}  // namespace seurat
