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

#ifndef VR_SEURAT_IMAGE_DEEP_EXR_IMAGE_READER_H_
#define VR_SEURAT_IMAGE_DEEP_EXR_IMAGE_READER_H_

#include <memory>
#include <vector>

#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "seurat/base/status.h"
#include "seurat/image/image_reader.h"

namespace seurat {
namespace image {

// This class implements deep EXR reading.
//
// Arbitrary strings are supported as channel names. The following reserved
// names are also accepted:
// * CONSTANT_ZERO: this channel is read with the constant value 0.0f.
// * CONSTANT_ONE: this channel is read with the constant value 1.0f.
class DeepExrImageReader : public ImageReader {
 public:
  // Create a DeepExrImageReader that reads from the deep EXR file contents in
  // |exr_contents|.  |exr_contents| must remain valid for the lifetime of this
  // DeepExrImageReader.
  static base::Status Create(absl::string_view exr_contents,
                             std::unique_ptr<DeepExrImageReader>* reader);
  ~DeepExrImageReader() override;

  // ImageReader implementation.
  int GetImageWidth() const override;
  int GetImageHeight() const override;
  void AddChannel(const Channel& channel) override;
  base::Status Read() const override;

  // Get the total count of samples in this deep EXR image.
  int GetSampleCount() const;

  // Get the per-pixel sample counts of the deep EXR image, in row-major order.
  absl::Span<int> GetSampleCounts() const;

 private:
  DeepExrImageReader();

  // Initializes this DeepExrImageReader with deep EXR file contents.
  base::Status Init(absl::string_view exr_contents);

  // DeepExrImageReader uses the PIMPL idiom to hide the internals of OpenEXR
  // from the public header.
  struct Impl;

  std::unique_ptr<Impl> impl_;
};

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_DEEP_EXR_IMAGE_READER_H_
