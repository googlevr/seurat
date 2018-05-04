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

#ifndef VR_SEURAT_IMAGE_EXR_UTIL_H_
#define VR_SEURAT_IMAGE_EXR_UTIL_H_

#include "absl/strings/string_view.h"
#include "OpenEXR/IlmImf/ImfGenericInputFile.h"
#include "OpenEXR/IlmImf/ImfIO.h"
#include "seurat/base/status.h"

namespace seurat {
namespace image {

// Adapt Imf IO to read from string objects, to match Google file
// {Get,Set}Contents APIs. OpenEXR only provides wrappers for Standard Library
// streams.
class ExrIStream : public Imf::IStream {
 public:
  explicit ExrIStream(absl::string_view bytestring)
      : Imf::IStream("BytestringIStream"),
        bytestring_(bytestring),
        file_position_(0) {}
  bool read(char c[], int n) override;
  Imath::Int64 tellg() override { return file_position_; }
  void seekg(Imath::Int64 pos) override { file_position_ = pos; }

 private:
  const absl::string_view bytestring_;  // The "file" contents.
  size_t file_position_;                // The "file" pointer.
};

// Adapt Imf IO to write to string objects, to match Google file
// {Get,Set}Contents APIs. OpenEXR only provides wrappers for Standard Library
// streams.
class ExrOStream : public Imf::OStream {
 public:
  explicit ExrOStream(std::string* bytestring)
      : Imf::OStream("BytestringOStream"),
        bytestring_(bytestring),
        file_position_(0) {}
  void write(const char c[], int n) override;
  Imath::Int64 tellp() override { return file_position_; }
  void seekp(Imath::Int64 pos) override;

 private:
  std::string* bytestring_;  // The "file" contents.
  size_t file_position_;     // The "file" pointer.
};

// Imf::DeepScanLineInputFile is the only Imf file type that cannot be read from
// a stream, without passing the header and version into the constructor. That's
// why we need this helper class to first read the header and version from the
// stream.
class ExrHeaderReader : Imf::GenericInputFile {
 public:
  // Reads the magic number, version and header from the |stream|.
  explicit ExrHeaderReader(Imf::IStream* stream) : stream_(stream) {
    readMagicNumberAndVersionField(*stream_, version_);
    header_.readFrom(*stream_, version_);
  }
  // Returns the version read from the stream (NOT the version of the header).
  int GetVersion() const { return version_; }
  // Returns the header.
  const Imf::Header& GetHeader() const { return header_; }

 private:
  Imf::IStream* stream_;
  int version_;
  Imf::Header header_;
};

// Return the image size from the EXR header.
base::Status ImageSizeFromHeader(const Imf::Header& header, int* image_width,
                                 int* image_height);

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_EXR_UTIL_H_
