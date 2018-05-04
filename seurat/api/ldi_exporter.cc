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

#include <memory>

#include "ion/math/vector.h"
#include "seurat/base/color.h"
#include "seurat/image/deep_exr_image_writer.h"

namespace seurat {
namespace api {

base::Status WriteDeepExr(int width, int height,
                          const std::vector<int>& sample_counts,
                          const std::vector<float>& rgba_colors,
                          const std::vector<float>& depths,
                          std::string* file_contents) {
  CHECK_EQ(width * height, sample_counts.size());
  CHECK_EQ(depths.size() * 4, rgba_colors.size());
  std::unique_ptr<image::DeepExrImageWriter> writer;
  SEURAT_RETURN_IF_ERROR(
      image::DeepExrImageWriter::Create(width, height, sample_counts, &writer));
  writer->AddChannel(
      {"R", &rgba_colors.data()[0], 4, image::ImageWriter::ValueType::kHalf});
  writer->AddChannel(
      {"G", &rgba_colors.data()[1], 4, image::ImageWriter::ValueType::kHalf});
  writer->AddChannel(
      {"B", &rgba_colors.data()[2], 4, image::ImageWriter::ValueType::kHalf});
  writer->AddChannel(
      {"A", &rgba_colors.data()[3], 4, image::ImageWriter::ValueType::kHalf});
  writer->AddChannel(
      {"Z", &depths.data()[0], 1, image::ImageWriter::ValueType::kFloat});
  return writer->Write(file_contents);
}

}  // namespace api
}  // namespace seurat
