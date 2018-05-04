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

#include "seurat/image/deep_exr_image_reader.h"

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
#include "OpenEXR/IlmImf/ImfDeepFrameBuffer.h"
#include "OpenEXR/IlmImf/ImfDeepScanLineInputFile.h"
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

// Internal state for DeepExrImageReader.
struct DeepExrImageReader::Impl {
  Impl() : image_width(0), image_height(0), sample_count(0) {}

  std::unique_ptr<ExrIStream> exr_istream;
  std::unique_ptr<Imf::DeepScanLineInputFile> file;
  std::unique_ptr<Imf::Header> header;
  int image_width;
  int image_height;
  std::vector<int> sample_counts;
  int sample_count;
  std::vector<Channel> channels;
};

// static
base::Status DeepExrImageReader::Create(
    absl::string_view exr_contents,
    std::unique_ptr<DeepExrImageReader>* reader) {
  std::unique_ptr<DeepExrImageReader> ret(new DeepExrImageReader());
  SEURAT_RETURN_IF_ERROR(ret->Init(exr_contents));
  *reader = std::move(ret);
  return base::OkStatus();
}

DeepExrImageReader::~DeepExrImageReader() {}

int DeepExrImageReader::GetImageWidth() const { return impl_->image_width; }

int DeepExrImageReader::GetImageHeight() const { return impl_->image_height; }

int DeepExrImageReader::GetSampleCount() const { return impl_->sample_count; }

absl::Span<int> DeepExrImageReader::GetSampleCounts() const {
  return absl::Span<int>(impl_->sample_counts.data(),
                         impl_->sample_counts.size());
}

void DeepExrImageReader::AddChannel(const Channel& channel) {
  impl_->channels.push_back(channel);
}

base::Status DeepExrImageReader::Read() const {
  if (!impl_) {
    return base::FailedPreconditionError("DeepExrImageReader::Init() failed");
  }

  // OpenEXR expects data for each channel of a deep framebuffer as a matrix of
  // pointers, laid out in row major order. The pointers for each pixel point to
  // an array of samples.  The sample count slice is unconditionally required.
  Imf::DeepFrameBuffer framebuffer;
  std::vector<unsigned int> sample_counts(impl_->image_width *
                                          impl_->image_height);
  framebuffer.insertSampleCountSlice(Imf::Slice(
      Imf::UINT, reinterpret_cast<char*>(sample_counts.data()),
      sizeof(decltype(sample_counts)::value_type),
      sizeof(decltype(sample_counts)::value_type) * impl_->image_width));

  // Iterate over all the channels and add a slice for each.
  const Imf::ChannelList& channels = impl_->header->channels();
  std::vector<std::vector<float*>> channel_ptr_list;
  for (const auto& channel : impl_->channels) {
    const absl::string_view& channel_name = channel.name;
    const bool channel_in_exr =
        (channels.findChannel(std::string(channel_name)) != nullptr);
    const bool is_reserved = IsReservedChannelName(channel_name);
    if (!channel_in_exr && !is_reserved) {
      return base::InvalidArgumentError(
          absl::Substitute("EXR channel \"$0\" missing", channel_name));
    } else if (channel_in_exr && is_reserved) {
      return base::InvalidArgumentError(absl::Substitute(
          "EXR reserved channel name \"$0\" cannot be used", channel_name));
    }

    std::vector<float*> channel_ptr;
    channel_ptr.reserve(impl_->sample_count);
    int sample_index = 0;
    for (auto sample_count : impl_->sample_counts) {
      int offset = sample_index * channel.element_period;
      channel_ptr.push_back(sample_count == 0 ? nullptr
                                              : channel.data + offset);
      sample_index += sample_count;
    }

    // Insert the slice for the channel.
    const int x_sampling = 1;
    const int y_sampling = 1;
    framebuffer.insert(
        std::string(channel_name),
        Imf::DeepSlice(Imf::PixelType::FLOAT,
                       reinterpret_cast<char*>(channel_ptr.data()),
                       sizeof(void*), sizeof(void*) * impl_->image_width,
                       channel.element_period * sizeof(float), x_sampling,
                       y_sampling, GetChannelDefaultValue(channel_name)));
    channel_ptr_list.emplace_back(std::move(channel_ptr));
  }

  try {
    // Read the sample counts and pixels for scanlines 0 through height-1.
    impl_->file->setFrameBuffer(framebuffer);
    impl_->file->readPixelSampleCounts(0, impl_->image_height - 1);
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

DeepExrImageReader::DeepExrImageReader() {}

base::Status DeepExrImageReader::Init(absl::string_view exr_contents) {
  auto impl = absl::make_unique<Impl>();

  impl->exr_istream = absl::make_unique<ExrIStream>(exr_contents);
  std::unique_ptr<ExrHeaderReader> header_reader;
  try {
    header_reader = absl::make_unique<ExrHeaderReader>(impl->exr_istream.get());
  } catch (const IEX_NAMESPACE::InputExc& input_exception) {
    return base::InvalidArgumentError(std::string("EXR ") +
                                      input_exception.what() +
                                      input_exception.stackTrace());
  }
  const int version = header_reader->GetVersion();
  impl->header = absl::make_unique<Imf::Header>(header_reader->GetHeader());
  impl->file = absl::make_unique<Imf::DeepScanLineInputFile>(
      *impl->header, impl->exr_istream.get(), version);

  // Validates the header and fetch the dimensions of the LDI in pixels.
  SEURAT_RETURN_IF_ERROR(ImageSizeFromHeader(
      *impl->header, &impl->image_width, &impl->image_height));

  // Setup a DeepFrameBuffer instance to read the sample counts.
  Imf::DeepFrameBuffer framebuffer;
  std::vector<unsigned int> sample_counts(impl->image_width *
                                          impl->image_height);
  framebuffer.insertSampleCountSlice(Imf::Slice(
      Imf::UINT, reinterpret_cast<char*>(sample_counts.data()),
      sizeof(decltype(sample_counts)::value_type),
      sizeof(decltype(sample_counts)::value_type) * impl->image_width));

  try {
    // Read the sample counts for scanlines 0 through height-1.
    impl->file->setFrameBuffer(framebuffer);
    impl->file->readPixelSampleCounts(0, impl->image_height - 1);
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

  // Count total number of samples.
  impl->sample_counts.reserve(sample_counts.size());
  for (const auto count : sample_counts) {
    impl->sample_counts.push_back(count);
    impl->sample_count += count;
  }

  // Initialization succeeded, so we now populate |impl_|.
  impl_ = std::move(impl);
  return base::OkStatus();
}

}  // namespace image
}  // namespace seurat
