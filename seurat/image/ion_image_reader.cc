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

#include "seurat/image/ion_image_reader.h"

#include "ion/image/conversionutils.h"
#include "ion/portgfx/glheaders.h"
#include "absl/strings/substitute.h"

namespace seurat {
namespace image {

namespace {

bool IsChannelNameSupported(absl::string_view channel_name) {
  std::array<std::string, 7> supported_names{
      {"R", "G", "B", "A", "CONSTANT_ZERO", "CONSTANT_ONE"}};
  if (std::find(supported_names.begin(), supported_names.end(), channel_name) !=
      supported_names.end()) {
    return true;
  } else {
    return false;
  }
}

// Checks a channel name to see if it is a reserved constant channel
bool IsConstantChannelName(absl::string_view name) {
  if (name == "CONSTANT_ZERO") return true;
  if (name == "CONSTANT_ONE") return true;
  return false;
}

// Returns the fill value for a reserved constant channel name, or 0.0f for all
// other channel names.
float GetChannelDefaultValue(absl::string_view name) {
  if (name == "CONSTANT_ZERO") return 0.0f;
  if (name == "CONSTANT_ONE") return 1.0f;
  return 0.0f;
}
// Checks a channel name to see if it is a reserved standard channel.
bool IsStandardChannelName(absl::string_view name) {
  if (name == "R") return true;
  if (name == "G") return true;
  if (name == "B") return true;
  if (name == "A") return true;
  return false;
}

// Returns the offset into an UInt8 pixel for a standard channel.
int GetChannelOffset(absl::string_view name) {
  if (name == "R") return 0;
  if (name == "G") return 1;
  if (name == "B") return 2;
  if (name == "A") return 3;
  LOG(ERROR) << name << " is not a standard channel name.";
  return 0;
}

}  // namespace

// static
base::Status IonImageReader::Create(absl::string_view image_contents,
                                    std::unique_ptr<IonImageReader>* reader) {
  std::unique_ptr<IonImageReader> ret(new IonImageReader());
  SEURAT_RETURN_IF_ERROR(ret->Init(image_contents));
  *reader = std::move(ret);
  return base::OkStatus();
}

int IonImageReader::GetImageWidth() const { return ion_image_->GetWidth(); }

int IonImageReader::GetImageHeight() const { return ion_image_->GetHeight(); }

void IonImageReader::AddChannel(const Channel& channel) {
  channels_.push_back(channel);
}

base::Status IonImageReader::Init(absl::string_view image_contents) {
  // Don't flip. Image flipping is handled by ingest.
  const bool kIonFlipImageVertically = false;
  // Don't wipe memory until the image is destroyed.
  const bool kIonDontWipe = false;

  // Decode into Ion Image.
  ion_image_ = ion::image::ConvertFromExternalImageData(
      image_contents.data(), image_contents.size(), kIonFlipImageVertically,
      kIonDontWipe, ion::base::AllocatorPtr());
  if (ion_image_.Get() == nullptr) {
    return base::InternalError("Decoding of image data failed.");
  }

  ion::gfx::Image::PixelFormat pixel_format =
      ion::gfx::Image::GetPixelFormat(ion_image_->GetFormat());
  switch (pixel_format.format) {
    case GL_RGB:
      num_channels_in_image_ = 3;
      break;
    case GL_RGBA:
      num_channels_in_image_ = 4;
      break;
    default:
      return base::InvalidArgumentError(
          "Pixel format must be either GL_RGBA or GL_RGB.");
  }

  if (pixel_format.type != GL_UNSIGNED_BYTE) {
    return base::InvalidArgumentError("Pixel type must be GL_UNSIGNED_BYTE.");
  }

  return base::OkStatus();
}

base::Status IonImageReader::Read() const {
  const int pixel_count = GetImageWidth() * GetImageHeight();
  const uint8* pixels = ion_image_->GetData()->GetData<uint8>();

  for (const Channel& channel : channels_) {
    if (!IsChannelNameSupported(channel.name)) {
      return base::InvalidArgumentError(
          absl::StrCat("Invalid channel name: ", channel.name));
    }
    if (channel.name == "A" && num_channels_in_image_ != 4) {
      return base::InvalidArgumentError("Image has no alpha channel.");
    }

    const uint8* pixel = pixels;
    float* element = channel.data;

    if (IsConstantChannelName(channel.name)) {
      // "CONSTANT_ZERO", "CONSTANT_ONE"
      const float value = GetChannelDefaultValue(channel.name);
      for (int i = 0; i < pixel_count; ++i) {
        *element = value;
        element += channel.element_period;
      }
    } else if (IsStandardChannelName(channel.name)) {
      // "R", "G", "B", "A"
      const int offset = GetChannelOffset(channel.name);
      DCHECK_GE(offset, 0);
      DCHECK_LE(offset, 3);
      for (int i = 0; i < pixel_count; ++i) {
        *element = pixel[offset] / 255.0f;
        element += channel.element_period;
        pixel += num_channels_in_image_;
      }
    }
  }

  return base::OkStatus();
}

}  // namespace image
}  // namespace seurat
