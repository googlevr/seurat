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

#include "seurat/base/structured_io.h"

#include <cstdint>

#include "ion/base/datacontainer.h"
#include "ion/base/logging.h"
#include "seurat/base/ion_util_no_gl.h"

namespace seurat {
namespace base {

using ion::math::Matrix4f;
using ion::math::Vector2i;

void StructureSink::WriteBytes(const char* bytes, size_t count) {
  sink_->Append(bytes, count);
}

void StructureSink::WriteString(const std::string& str) {
  WritePod<uint64_t>(str.size());
  WriteBytes(str.c_str(), str.size());
}

void StructureSink::WriteImage(ion::gfx::ImagePtr image) {
  CHECK(image.Get());
  WritePod<int>(image->GetWidth());
  WritePod<int>(image->GetHeight());
  WritePod<ion::gfx::Image::Format>(image->GetFormat());
  WriteBytes(image->GetData()->GetData<char>(), image->GetDataSize());
}

void StructureSource::ReadBytes(char* bytes_out, size_t count) {
  CHECK_GE(source_->Available(), count);
  UncheckedArrayByteSink outsink(bytes_out);
  source_->CopyTo(&outsink, count);
}

std::string StructureSource::ReadString() {
  uint64_t length = ReadPod<uint64_t>();
  CHECK_GE(source_->Available(), length);
  std::string str;
  StringByteSink strsink(&str);
  source_->CopyTo(&strsink, length);
  return std::string(str);
}

ion::gfx::ImagePtr StructureSource::ReadImage() {
  const int width = ReadPod<int>();
  const int height = ReadPod<int>();
  const auto format = ReadPod<ion::gfx::Image::Format>();
  ion::gfx::ImagePtr image =
      seurat::base::CreateImage(format, {width, height});
  ion::base::DataContainerPtr container = image->GetData();
  ReadBytes(container->GetMutableData<char>(), image->GetDataSize());
  return image;
}

}  // namespace base
}  // namespace seurat
