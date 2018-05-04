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

#include "seurat/artifact/exr_exporter.h"

#include "absl/strings/str_cat.h"
#include "seurat/base/array2d_util.h"
#include "seurat/base/status.h"
#include "seurat/image/exr_image_writer.h"

namespace seurat {
namespace artifact {

using image::ExrImageWriter;
using image::Image4f;

base::Status ExrExporter::Process(Artifact* artifact) const {
  if (!artifact->texture) {
    return base::UnimplementedError("No texture to export.");
  }
  Image4f image = *artifact->texture;

  std::string filename = absl::StrCat(basename_, ".exr");
  // Flip along the Y-axis.
  base::FlipArrayVertical(&image);
  // Write the Image to a string.
  std::string encoded;

  const auto color_type = ExrImageWriter::ValueType::kFloat;
  std::unique_ptr<ExrImageWriter> writer;
  SEURAT_RETURN_IF_ERROR(
      ExrImageWriter::Create(image.Width(), image.Height(), &writer));
  const int dimension = Image4f::ElementType::kDimension;
  writer->AddChannel({"R", &image.data()[0][0], dimension, color_type});
  writer->AddChannel({"G", &image.data()[0][1], dimension, color_type});
  writer->AddChannel({"B", &image.data()[0][2], dimension, color_type});
  writer->AddChannel({"A", &image.data()[0][3], dimension, color_type});
  SEURAT_RETURN_IF_ERROR(writer->Write(&encoded));
  return filesystem_->SetContents(filename, encoded);
}

}  // namespace artifact
}  // namespace seurat
