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

#ifndef VR_SEURAT_COMPONENT_ICE_IO_H_
#define VR_SEURAT_COMPONENT_ICE_IO_H_

#include <memory>

#include "seurat/base/structured_io.h"
#include "seurat/component/component.h"

namespace seurat {
namespace component {

// Deserializes an ICE file from the |stream| into a Component hierarchy. The
// stream must start with an IceHeader struct, e.g. as written by WriteIce.
// Returns the deserialized Component hierarchy, or nullptr in case of an
// invalid header or unrecognized format.
std::unique_ptr<const Component> ReadIce(base::StructureSource* stream);

// Serializes a Component hierarchy into an ICE file format |stream|. The stream
// will start with an IceHeader struct.
void WriteIce(const Component& source, base::StructureSink* stream);

// ICE Format API Implementation details live here.
namespace ice_io_internal {

// Defines ICE file serialized data format version numbers. Higher numbers are
// newer versions and lower are older.
//
// WARNING! Implement this number every time the serialized format for any
// Component changes.
static constexpr int kCurrent = 1;

// Defines ICE Format parsing for API version |Revision|. This in particular
// allows testing different versions of the API against each other.
class IceFormat {
 public:
  explicit IceFormat(int revision) : revision_(revision) {}

  // Deserializes and validates the ICE format  header from the |stream|.
  // Returns true if the header encodes a valid ICE signature. The stream must
  // start with an IceHeader struct, e.g. as written by WriteIce.
  bool ReadIceHeader(base::StructureSource* stream);

  // Stores a current-version ICE format header to the |stream|.
  void WriteIceHeader(base::StructureSink* stream);

  // Deserializes an ICE file from the |stream| into a Component hierarchy. The
  // stream must start with an IceHeader struct, e.g. as written by WriteIce.
  // Returns the deserialized Component hierarchy, or nullptr in case of an
  // invalid header or unrecognized format.
  std::unique_ptr<const Component> ReadIce(base::StructureSource* stream);

  // Serializes a Component hierarchy into an ICE file format |stream|. The
  // stream will start with an IceHeader struct.
  void WriteIce(const Component& source, base::StructureSink* stream);

 private:
  // Defines a POD header for an ICE file.
  struct IceFormatHeader {
    // This member must contain 'IceF' to indicate the Ice Format.
    uint8_t signature_fourcc[4];
    // Stores the version of the format that serialized the data.
    uint32_t format_version;
  };

  // Holds the version number of this instance of the format API. See
  // ice_io_internal::kCurrent.
  const int revision_;

  enum class HeaderStatus {
    kCorrect,
    kInvalidSignature,
    kObsoleteVersion,
    kNewerUnknownVersion,
  };
  HeaderStatus IsValidHeader(const IceFormatHeader& header) const;
};

}  // namespace ice_io_internal

}  // namespace component
}  // namespace seurat

#endif  // VR_SEURAT_COMPONENT_ICE_IO_H_
