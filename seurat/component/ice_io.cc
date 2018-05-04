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

#include "seurat/component/ice_io.h"

#include <array>

namespace seurat {
namespace component {

using ice_io_internal::IceFormat;

std::unique_ptr<const Component> ReadIce(base::StructureSource* stream) {
  IceFormat ice_format(ice_io_internal::kCurrent);
  return ice_format.ReadIce(stream);
}

void WriteIce(const Component& source, base::StructureSink* stream) {
  IceFormat ice_format(ice_io_internal::kCurrent);
  ice_format.WriteIce(source, stream);
}

namespace ice_io_internal {

std::unique_ptr<const Component> IceFormat::ReadIce(
    base::StructureSource* stream) {
  IceFormatHeader header;
  header = stream->ReadPod<IceFormatHeader>();
  std::unique_ptr<const Component> result;
  HeaderStatus status = IsValidHeader(header);
  if (status == HeaderStatus::kCorrect) {
    result = Component::Create(stream);
  } else {
    // TODO(ernstm) Return this error message inside the base::Status object
    // once we make that conversion.
    static std::array<const char*, 4> status_messages{
        {nullptr, "bad signature", "old version", " unrecognized new version"}};
    LOG(ERROR) << "Invalid ICE header: "
               << status_messages[static_cast<int>(status)];
  }
  return result;
}

void IceFormat::WriteIce(const Component& source, base::StructureSink* stream) {
  IceFormatHeader header{{'I', 'c', 'e', 'F'},
                         static_cast<uint32_t>(revision_)};
  stream->WritePod(header);
  source.Write(stream);
}

IceFormat::HeaderStatus IceFormat::IsValidHeader(
    const IceFormatHeader& header) const {
  HeaderStatus status = HeaderStatus::kCorrect;
  if (!(header.signature_fourcc[0] == 'I' &&
        header.signature_fourcc[1] == 'c' &&
        header.signature_fourcc[2] == 'e' &&
        header.signature_fourcc[3] == 'F')) {
    status = HeaderStatus::kInvalidSignature;
  } else if (header.format_version < static_cast<uint32_t>(revision_)) {
    status = HeaderStatus::kObsoleteVersion;
  } else if (header.format_version > static_cast<uint32_t>(revision_)) {
    status = HeaderStatus::kNewerUnknownVersion;
  }
  return status;
}

}  // namespace ice_io_internal

}  // namespace component
}  // namespace seurat
