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

#ifndef VR_SEURAT_IMAGE_IMAGE_READER_H_
#define VR_SEURAT_IMAGE_IMAGE_READER_H_

#include "absl/strings/string_view.h"
#include "seurat/base/status.h"

namespace seurat {
namespace image {

// Instances of implementations of ImageReader are created through a factory
// function, and Channel definitions are incrementally added to it with
// AddChannel().  The image data is then read into the channels with Read().
class ImageReader {
 public:
  struct Channel {
    // The name of this channel, to be read from the image file. Some
    // implementations may support reserved names for constant values or indexed
    // channels.
    absl::string_view name;

    // The output buffer which the data in the image file will be read into.
    float* data;

    // The period of the elements in this channel.  If '1', the elements are
    // written contiguously into |data|; if |2|, the elements are written to
    // every other index of |data|, etc.
    int element_period;
  };

  virtual ~ImageReader() = default;

  // Returns the width of the image.
  virtual int GetImageWidth() const = 0;

  // Returns the height of the image.
  virtual int GetImageHeight() const = 0;

  // Add a channel definition.
  //
  // TODO(ernstm): Should we add all channels with Create() and then validate
  // them? Probably yes.
  virtual void AddChannel(const Channel& channel) = 0;

  // Read the image into the defined channels.
  virtual base::Status Read() const = 0;
};

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_IMAGE_READER_H_
