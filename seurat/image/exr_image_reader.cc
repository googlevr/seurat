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

#include "seurat/image/exr_image_reader.h"

#include <algorithm>
#include <limits>
#include <string>

#include "absl/memory/memory.h"
#include "absl/strings/string_view.h"
#include "absl/strings/substitute.h"
#include "IlmBase/Half/half.h"
#include "IlmBase/Iex/IexBaseExc.h"
#include "IlmBase/Iex/IexMacros.h"
#include "OpenEXR/IlmImf/ImfChannelList.h"
#include "OpenEXR/IlmImf/ImfIO.h"
#include "OpenEXR/IlmImf/ImfInputFile.h"
#include "OpenEXR/IlmImf/ImfInt64.h"
#include "OpenEXR/IlmImf/ImfPartType.h"
#include "seurat/image/exr_util.h"

namespace seurat {
namespace image {
namespace {

// Checks a channel name to see if it is reserved.
bool IsReservedChannelName(absl::string_view name) {
  if (name == "CONSTANT_ZERO") {
    return true;
  }
  if (name == "CONSTANT_ONE") {
    return true;
  }
  return false;
}

// Returns the fill value for a "reserved" channel name.
float GetChannelDefaultValue(absl::string_view name) {
  if (name == "CONSTANT_ZERO") {
    return 0.0f;
  }
  if (name == "CONSTANT_ONE") {
    return 1.0f;
  }
  return 0.0f;
}

}  // namespace

// Internal state for ExrImageReader.
struct ExrImageReader::Impl {
  Impl() : image_width(0), image_height(0) {}

  std::unique_ptr<ExrIStream> exr_istream;
  std::unique_ptr<Imf::InputFile> file;
  std::unique_ptr<Imf::Header> header;
  int image_width;
  int image_height;
  std::vector<Channel> channels;
};

// static
base::Status ExrImageReader::Create(absl::string_view exr_contents,
                                    std::unique_ptr<ExrImageReader>* reader) {
  std::unique_ptr<ExrImageReader> ret(new ExrImageReader());
  SEURAT_RETURN_IF_ERROR(ret->Init(exr_contents));
  *reader = std::move(ret);
  return base::OkStatus();
}

ExrImageReader::~ExrImageReader() {}

int ExrImageReader::GetImageWidth() const { return impl_->image_width; }

int ExrImageReader::GetImageHeight() const { return impl_->image_height; }

void ExrImageReader::AddChannel(const Channel& channel) {
  impl_->channels.push_back(channel);
}

base::Status ExrImageReader::Read() const {
  if (!impl_) {
    return base::FailedPreconditionError("ExrImageReader::Init() failed");
  }

  Imf::FrameBuffer framebuffer;
  // Iterate over all the channels and add a slice for each.
  const Imf::ChannelList& channels = impl_->header->channels();
  std::vector<std::vector<float*>> channel_ptr_list;
  for (const auto& channel : impl_->channels) {
    const absl::string_view& channel_name = channel.name;
    const bool is_reserved = IsReservedChannelName(channel_name);
    bool channel_in_exr = false;

    // Note the OpenEXR library documents this can throw ArgExc, but the
    // implementation doesn't actually throw. We've wrapped it in case it's
    // fixed.
    try {
      channel_in_exr =
          (channels.findChannel(std::string(channel_name)) != nullptr);
    } catch (const IEX_NAMESPACE::ArgExc& arg_exception) {
      // This exception is not an error. It only indicates that no channel
      // with the given name was found.
      channel_in_exr = false;
    } catch (...) {
      return base::InternalError("Unrecognized EXR exception");
    }

    if (!channel_in_exr && !is_reserved) {
      std::string error_message = absl::Substitute(
          "EXR channel \"$0\" missing. Existing channels: ", channel_name);
      Imf::ChannelList::ConstIterator iter = channels.begin();
      while (iter != channels.end()) {
        error_message += iter.name();
        ++iter;
        if (iter != channels.end()) {
          error_message += std::string(", ");
        }
      }
      return base::InvalidArgumentError(error_message);
    } else if (channel_in_exr && is_reserved) {
      return base::InvalidArgumentError(absl::Substitute(
          "EXR reserved channel name \"$0\" cannot be used", channel_name));
    }

    // Insert the slice for the channel.
    const int x_sampling = 1;
    const int y_sampling = 1;
    const size_t x_stride = channel.element_period * sizeof(float);
    const size_t y_stride = x_stride * impl_->image_width;
    framebuffer.insert(
        std::string(channel_name),
        Imf::Slice(Imf::PixelType::FLOAT, reinterpret_cast<char*>(channel.data),
                   x_stride, y_stride, x_sampling, y_sampling,
                   GetChannelDefaultValue(channel_name)));
  }

  // TODO(b/36897469): Improve test coverage for error cases.
  try {
    // Read the pixels for scanlines 0 through height-1.
    impl_->file->setFrameBuffer(framebuffer);
    impl_->file->readPixels(0, impl_->image_height - 1);
  } catch (const IEX_NAMESPACE::ArgExc& arg_exception) {
    return base::InvalidArgumentError(std::string("EXR ") +
                                      arg_exception.what() +
                                      arg_exception.stackTrace());
  } catch (const IEX_NAMESPACE::IoExc& io_exception) {
    return base::OutOfRangeError(std::string("EXR ") + io_exception.what() +
                                 io_exception.stackTrace());
  } catch (...) {
    return base::InternalError("Unrecognized EXR exception");
  }
  return base::OkStatus();
}

ExrImageReader::ExrImageReader() {}

base::Status ExrImageReader::Init(absl::string_view exr_contents) {
  auto impl = absl::make_unique<Impl>();

  // TODO(b/36897469): Improve test coverage for error cases.
  try {
    impl->exr_istream = absl::make_unique<ExrIStream>(exr_contents);
    // Construct the InputFile from the ExrIStream.
    impl->file = absl::make_unique<Imf::InputFile>(*impl->exr_istream);
    impl->header = absl::make_unique<Imf::Header>(impl->file->header());
  } catch (const IEX_NAMESPACE::ArgExc& arg_exception) {
    return base::InvalidArgumentError(std::string("EXR ") +
                                      arg_exception.what() +
                                      arg_exception.stackTrace());
  } catch (const IEX_NAMESPACE::IoExc& io_exception) {
    return base::OutOfRangeError(std::string("EXR ") + io_exception.what() +
                                 io_exception.stackTrace());
  } catch (...) {
    return base::InternalError("Unrecognized EXR exception");
  }

  // Validates the header and fetch the dimensions of the LDI in pixels.
  SEURAT_RETURN_IF_ERROR(ImageSizeFromHeader(*impl->header, &impl->image_width,
                                             &impl->image_height));

  // Initialization succeeded, so we now populate |impl_|.
  impl_ = std::move(impl);
  return base::OkStatus();
}

}  // namespace image
}  // namespace seurat
