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

#include "seurat/ingest/ldi_loader.h"

#include <string>
#include <thread>  // NOLINT

#include "absl/strings/substitute.h"
#include "seurat/base/array2d_util.h"
#include "seurat/base/color.h"
#include "seurat/base/status_util.h"
#include "seurat/image/deep_exr_image_reader.h"
#include "seurat/image/exr_image_reader.h"
#include "seurat/image/image.h"
#include "seurat/image/ion_image_reader.h"

namespace seurat {
namespace ingest {
namespace {

base::Status LoadLdiFile(const api::proto::LdiFile& proto,
                         base::FileSystem* file_system, image::Ldi4f* ldi) {
  // If an error occurs, this extension will be added to provide additional
  // context. Construct it here once and use for all error cases.
  const base::Status status_extension = base::InvalidArgumentError(
      absl::StrCat("Failed loading Ldi file ", proto.path()));

  std::string exr_contents;
  SEURAT_RETURN_IF_ERROR_EXTENDED(
      file_system->GetContents(proto.path(), &exr_contents), status_extension);

  std::unique_ptr<image::DeepExrImageReader> reader;
  SEURAT_RETURN_IF_ERROR_EXTENDED(
      image::DeepExrImageReader::Create(exr_contents, &reader),
      status_extension);

  const ion::math::Vector2i size(reader->GetImageWidth(),
                                 reader->GetImageHeight());
  const int sample_count = reader->GetSampleCount();
  std::vector<int> sample_counts;
  // Reserve space for one additional element to avoid a re-allocation in the
  // Ldi constructor, when sample counts are converted to offsets with a
  // sentinel at the end.
  sample_counts.reserve(size[0] * size[1] + 1);
  std::copy(reader->GetSampleCounts().begin(), reader->GetSampleCounts().end(),
            std::back_inserter(sample_counts));
  std::vector<base::Color4f> colors(sample_count);
  std::vector<float> depths(sample_count);

  reader->AddChannel({proto.channel_0(), &colors.data()[0][0], 4});
  reader->AddChannel({proto.channel_1(), &colors.data()[0][1], 4});
  reader->AddChannel({proto.channel_2(), &colors.data()[0][2], 4});
  reader->AddChannel({proto.channel_alpha(), &colors.data()[0][3], 4});
  reader->AddChannel({proto.channel_depth(), &depths.data()[0], 1});
  SEURAT_RETURN_IF_ERROR_EXTENDED(reader->Read(), status_extension);

  *ldi = image::Ldi4f(size, std::move(sample_counts), std::move(colors),
                      std::move(depths));
  return base::OkStatus();
}

base::Status CreateImageReader(absl::string_view encoded_image,
                               std::unique_ptr<image::ImageReader>* reader) {
  // Try the ExrImageReader first. Use it if encoded_image is a valid Exr file.
  std::unique_ptr<image::ExrImageReader> exr_reader;
  base::Status exr_status =
      image::ExrImageReader::Create(encoded_image, &exr_reader);
  if (exr_status.ok()) {
    *reader = std::move(exr_reader);
    return exr_status;
  }

  // Now try the IonImageReader. Use it if Ion can reader the encoded_image.
  std::unique_ptr<image::IonImageReader> ion_reader;
  base::Status ion_status =
      image::IonImageReader::Create(encoded_image, &ion_reader);
  if (ion_status.ok()) {
    *reader = std::move(ion_reader);
    return ion_status;
  }

  // At this point, none of the image readers could handle the
  // encoded_image. Return the statuses of all readers and a high-level error
  // message.
  base::Status status;
  base::UpdateStatus(&status, exr_status);
  base::UpdateStatus(&status, ion_status);
  base::UpdateStatus(&status,
                     base::InvalidArgumentError(
                         "Image file not recognized by any image reader!"));
  return status;
}

base::Status LoadImage4File(const api::proto::Image4File& proto,
                            base::FileSystem* file_system,
                            image::Image4f* image) {
  // If an error occurs, this extension will be added to provide additional
  // context. Construct it here once and use for all error cases.
  const base::Status status_extension = base::InvalidArgumentError(
      absl::StrCat("Failed loading Image4 file ", proto.path()));

  std::string encoded_image;
  SEURAT_RETURN_IF_ERROR_EXTENDED(
      file_system->GetContents(proto.path(), &encoded_image), status_extension);

  std::unique_ptr<image::ImageReader> reader;
  SEURAT_RETURN_IF_ERROR_EXTENDED(CreateImageReader(encoded_image, &reader),
                                     status_extension);

  image->Resize(reader->GetImageWidth(), reader->GetImageHeight());
  const int dimension = image::Image4f::ElementType::kDimension;
  reader->AddChannel({proto.channel_0(), &image->data()[0][0], dimension});
  reader->AddChannel({proto.channel_1(), &image->data()[0][1], dimension});
  reader->AddChannel({proto.channel_2(), &image->data()[0][2], dimension});
  reader->AddChannel({proto.channel_alpha(), &image->data()[0][3], dimension});
  SEURAT_RETURN_IF_ERROR_EXTENDED(reader->Read(), status_extension);

  // Flip the image along the Y-axis.
  base::FlipArrayVertical(image);

  return base::OkStatus();
}

base::Status LoadImage1File(const api::proto::Image1File& proto,
                            base::FileSystem* file_system,
                            image::Image1f* image) {
  // If an error occurs, this extension will be added to provide additional
  // context. Construct it here once and use for all error cases.
  const base::Status status_extension = base::InvalidArgumentError(
      absl::StrCat("Failed loading Image1 file ", proto.path()));

  std::string encoded_image;
  SEURAT_RETURN_IF_ERROR_EXTENDED(
      file_system->GetContents(proto.path(), &encoded_image), status_extension);

  std::unique_ptr<image::ImageReader> reader;
  SEURAT_RETURN_IF_ERROR_EXTENDED(CreateImageReader(encoded_image, &reader),
                                     status_extension);
  image->Resize(reader->GetImageWidth(), reader->GetImageHeight());
  const int dimension = image::Image1f::ElementType::kDimension;
  reader->AddChannel({proto.channel_0(), &image->data()[0][0], dimension});
  SEURAT_RETURN_IF_ERROR_EXTENDED(reader->Read(), status_extension);

  // Flip the image along the Y-axis.
  base::FlipArrayVertical(image);

  return base::OkStatus();
}

base::Status LoadDepthImageFile(const api::proto::DepthImageFile& proto,
                                base::FileSystem* file_system,
                                image::Ldi4f* ldi) {
  image::Image4f color_image;
  image::Image1f depth_image;
  // Spin up another thread to load and decode the color and depth images in
  // parallel.
  base::Status image1_status;
  base::Status image4_status;
  std::thread image4_thread([&]() {
    image4_status = LoadImage4File(proto.color(), file_system, &color_image);
  });
  image1_status = LoadImage1File(proto.depth(), file_system, &depth_image);
  // WARNING! All code that can return from this function must be after the
  // join!
  image4_thread.join();

  // Handle errors after both loads are complete. Otherwise, we might return
  // before image4_thread has written its return value to image4_status,
  // resulting in a crash.
  SEURAT_RETURN_IF_ERROR(image1_status);
  SEURAT_RETURN_IF_ERROR(image4_status);

  if (color_image.GetSize() != depth_image.GetSize()) {
    return base::InvalidArgumentError(
        absl::Substitute("Color image and depth image differ in size. Color: "
                         "$0 x $1. Depth: $2 x $3.",
                         color_image.Width(), color_image.Height(),
                         depth_image.Width(), depth_image.Height()));
  }
  const ion::math::Vector2i size(color_image.GetSize());

  std::vector<int> sample_counts;
  std::vector<base::Color4f> colors;
  std::vector<float> depths;

  for (int y = 0; y < size[1]; ++y) {
    for (int x = 0; x < size[0]; ++x) {
      if (color_image.At(x, y)[3] != 0.0f) {
        colors.push_back(color_image.At(x, y));
        depths.push_back(depth_image.At(x, y)[0]);
        sample_counts.push_back(1);
      } else {
        sample_counts.push_back(0);
      }
    }
  }

  *ldi = image::Ldi4f(size, std::move(sample_counts), std::move(colors),
                      std::move(depths));
  return base::OkStatus();
}

}  // namespace

base::Status LdiLoader::Load(const api::proto::Ldi& proto,
                             base::FileSystem* file_system,
                             image::Ldi4f* ldi) const {
  switch (proto.ldi_case()) {
    case api::proto::Ldi::LdiCase::LDI_NOT_SET:
      return base::InvalidArgumentError("Ldi.ldi not set");
    case api::proto::Ldi::LdiCase::kLdiFile:
      return LoadLdiFile(proto.ldi_file(), file_system, ldi);
    case api::proto::Ldi::LdiCase::kDepthImageFile:
      return LoadDepthImageFile(proto.depth_image_file(), file_system, ldi);
  }
}

}  // namespace ingest
}  // namespace seurat
