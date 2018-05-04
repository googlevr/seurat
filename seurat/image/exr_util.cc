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

#include "seurat/image/exr_util.h"

#include "absl/strings/substitute.h"
#include "IlmBase/Iex/IexBaseExc.h"
#include "IlmBase/Iex/IexMacros.h"

namespace seurat {
namespace image {

bool ExrIStream::read(char c[], int n) {
  const size_t bytestring_size = bytestring_.size();
  const size_t available_bytes = bytestring_size - file_position_;
  // Use exceptions for error handling to follow the specification of
  // Imf::IStream::read.  ASSERT is an OpenEXR library macro that throws an
  // exception.
  ASSERT(c != nullptr, IEX_NAMESPACE::IoExc, "null input buffer");
  ASSERT(n >= 0, IEX_NAMESPACE::IoExc, "negative read count");
  ASSERT(file_position_ <= bytestring_size, IEX_NAMESPACE::IoExc,
         "file position past end");
  ASSERT(static_cast<size_t>(n) <= available_bytes, IEX_NAMESPACE::IoExc,
         "read count greater than available bytes");
  bytestring_.copy(c, n, file_position_);
  file_position_ += n;
  return file_position_ < bytestring_size;
}

void ExrOStream::write(const char c[], int n) {
  // Use exceptions for error handling to follow the specification of
  // Imf::OStream::write.  ASSERT is an OpenEXR library macro that throws an
  // exception.
  ASSERT(bytestring_ != nullptr, IEX_NAMESPACE::IoExc, "null output string");
  ASSERT(c != nullptr, IEX_NAMESPACE::IoExc, "null output buffer");
  ASSERT(n >= 0, IEX_NAMESPACE::IoExc, "negative write count");
  const size_t bytestring_size = bytestring_->size();
  ASSERT(file_position_ <= bytestring_size, IEX_NAMESPACE::IoExc,
         "file position past end");
  const size_t end_fp = file_position_ + n;
  if (end_fp > bytestring_->size()) {
    bytestring_->resize(end_fp);
  }
  bytestring_->replace(file_position_, n, c, n);
  file_position_ = end_fp;
}

void ExrOStream::seekp(Imath::Int64 pos) {
  ASSERT(bytestring_ != nullptr, IEX_NAMESPACE::IoExc, "null output string");
  file_position_ = pos;
  // ImfIO.h doesn't specify whether seeking past EOF is valid. We assume that
  // it is and it grows the output file accordingly.
  if (file_position_ > bytestring_->size()) {
    bytestring_->resize(file_position_);
  }
}

base::Status ImageSizeFromHeader(const Imf::Header& header, int* image_width,
                                 int* image_height) {
  // Read the size of the LDI from the EXR header.
  const Imath::Box2i& data_window = header.dataWindow();
  const Imath::Box2i& display_window = header.displayWindow();
  if (data_window.min.x != 0 || data_window.min.y != 0) {
    return base::InvalidArgumentError(
        absl::Substitute("EXR data_window.min=($0, $1) != (0, 0)",
                         data_window.min.x, data_window.min.y));
  }
  if (display_window.min.x != 0 || display_window.min.y != 0) {
    return base::InvalidArgumentError(
        absl::Substitute("EXR display_window.min=($0, $1) != (0, 0)",
                         display_window.min.x, display_window.min.y));
  }
  if (data_window.max != display_window.max) {
    return base::InvalidArgumentError(absl::Substitute(
        "EXR data_window.max=($0, $1) != display_window.max=($2, $3)",
        data_window.max.x, data_window.max.y, display_window.max.x,
        display_window.max.y));
  }

  *image_width = data_window.max.x - data_window.min.x + 1;
  *image_height = data_window.max.y - data_window.min.y + 1;
  return base::OkStatus();
}

}  // namespace image
}  // namespace seurat
