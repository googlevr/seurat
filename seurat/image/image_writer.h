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

#ifndef VR_SEURAT_IMAGE_IMAGE_WRITER_H_
#define VR_SEURAT_IMAGE_IMAGE_WRITER_H_

#include <memory>
#include <string>
#include <vector>

#include "absl/strings/string_view.h"
#include "seurat/base/status.h"

namespace seurat {
namespace image {

// Serializes selected channels of an image into a file format defined by each
// implementation, storing the resulting data blob in a string.
// Implementation-provided factory function Create constructs instances
// of the writer. Channel definitions are incrementally added to it with
// AddChannel(). The image data is then written using Write().
class ImageWriter {
 public:
  // An enum for each value type supported. Some of these will result
  // in loss of precision and not all of them are supported by all
  // ImageWriter implementations.
  enum ValueType {
    kUnspecified = 0,
    kHalf = 1,
    kFloat = 2,
  };

  // Represents one "channel" of data.  All the pixels (or samples for deep
  // images) are written contiguously from a data buffer.  Each Channel
  // represent how one data component of the pixel (or sample) is to be arranged
  // in the input.
  struct Channel {
    // The name of this channel, to be written to the file.
    absl::string_view name;

    // The input buffer which will be written into the file.
    const float* data;

    // The period of the elements in this channel.  If '1', the elements are
    // read contiguously from |data|; if |2|, the elements are read from every
    // other index of |data|, etc.  If |0|, the first index in |data| is read
    // for every element.
    int element_period;

    // The desired output type of this channel in the file.
    ValueType out_type;
  };

  virtual ~ImageWriter() = default;

  // Add a channel definition.
  virtual void AddChannel(const Channel& channel) = 0;

  // Write the image.
  virtual base::Status Write(std::string* exr_contents) const = 0;
};

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_IMAGE_WRITER_H_
