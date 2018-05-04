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

#ifndef VR_SEURAT_IMAGE_ION_IMAGE_READER_H_
#define VR_SEURAT_IMAGE_ION_IMAGE_READER_H_

#include "ion/gfx/image.h"
#include "absl/strings/string_view.h"
#include "seurat/base/status.h"
#include "seurat/image/image_reader.h"

namespace seurat {
namespace image {

// An image reader implementation that uses Ion's image I/O internally.
//
// Only the following reserved channel names are accepted:
// * R: The channel with index 0.
// * G: The channel with index 1.
// * B: The channel with index 2.
// * A: The channel with index 3.
// * CONSTANT_ZERO: This channel is read with the constant value 0.0f.
// * CONSTANT_ONE: This channel is read with the constant value 1.0f.
class IonImageReader : public ImageReader {
 public:
  // Create an IonImageReader that reads from the file contents in
  // |image_contents|.  |image_contents| must remain valid for the lifetime of
  // this IonImageReader.
  static base::Status Create(absl::string_view image_contents,
                             std::unique_ptr<IonImageReader>* reader);
  ~IonImageReader() override = default;

  // ImageReader implementation.
  int GetImageWidth() const override;
  int GetImageHeight() const override;
  void AddChannel(const Channel& channel) override;
  base::Status Read() const override;

 private:
  IonImageReader() = default;

  // Initializes this ImageReader with image file contents.
  base::Status Init(absl::string_view image_contents);

  // The Ion image read from the file during initialization.
  ion::gfx::ImagePtr ion_image_;

  // Number of channels in the Ion image.
  int num_channels_in_image_;

  // The list of channels to be read.
  std::vector<Channel> channels_;
};

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_ION_IMAGE_READER_H_
