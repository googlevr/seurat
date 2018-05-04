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

#include "seurat/image/deep_exr_image_writer.h"

#include <numeric>

#include "absl/strings/string_view.h"
#include "absl/strings/substitute.h"
#include "IlmBase/Half/half.h"
#include "IlmBase/Iex/IexBaseExc.h"
#include "IlmBase/Iex/IexMacros.h"
#include "OpenEXR/IlmImf/ImfChannelList.h"
#include "OpenEXR/IlmImf/ImfDeepFrameBuffer.h"
#include "OpenEXR/IlmImf/ImfDeepScanLineOutputFile.h"
#include "OpenEXR/IlmImf/ImfIO.h"
#include "OpenEXR/IlmImf/ImfInt64.h"
#include "OpenEXR/IlmImf/ImfOutputFile.h"
#include "OpenEXR/IlmImf/ImfPartType.h"
#include "seurat/image/exr_util.h"

namespace seurat {
namespace image {

// static
base::Status DeepExrImageWriter::Create(
    int image_width, int image_height, std::vector<int> sample_counts,
    std::unique_ptr<DeepExrImageWriter>* writer) {
  std::unique_ptr<DeepExrImageWriter> ret(new DeepExrImageWriter());
  SEURAT_RETURN_IF_ERROR(
      ret->Init(image_width, image_height, std::move(sample_counts)));
  *writer = std::move(ret);
  return base::OkStatus();
}

DeepExrImageWriter::~DeepExrImageWriter() {}

void DeepExrImageWriter::AddChannel(const Channel& channel) {
  channels_.push_back(channel);
}

base::Status DeepExrImageWriter::Write(std::string* exr_contents) const {
  // Creates an EXR header from the given channels.
  const float kAspectRatio = 1.0f;
  const Imath::V2f kWindowCenter(0.0f, 0.0f);
  const float kWindowWidth = 1.0f;
  Imf::Header header(image_width_, image_height_, kAspectRatio, kWindowCenter,
                     kWindowWidth, Imf::INCREASING_Y, Imf::ZIPS_COMPRESSION);
  header.setType(Imf::DEEPSCANLINE);

  const int pixel_count = image_width_ * image_height_;
  const int sample_count =
      std::accumulate(sample_counts_.begin(), sample_counts_.end(), 0);

  // Create the DeepFramebuffer and first insert the sample count slice.
  Imf::DeepFrameBuffer framebuffer;
  std::vector<unsigned int> sample_counts(sample_counts_.begin(),
                                          sample_counts_.end());
  framebuffer.insertSampleCountSlice(
      Imf::Slice(Imf::UINT, reinterpret_cast<char*>(sample_counts.data()),
                 sizeof(decltype(sample_counts)::value_type),
                 sizeof(decltype(sample_counts)::value_type) * image_width_));

  // Iterate over all the channels and add a slice for each.
  std::vector<std::vector<half>> half_channel_list;
  std::vector<std::vector<const half*>> half_channel_ptr_list;
  std::vector<std::vector<const float*>> float_channel_ptr_list;
  for (const Channel& channel : channels_) {
    Imf::PixelType pixel_type = Imf::PixelType::NUM_PIXELTYPES;
    switch (channel.out_type) {
      case ValueType::kHalf:
        pixel_type = Imf::PixelType::HALF;
        break;
      case ValueType::kFloat:
        pixel_type = Imf::PixelType::FLOAT;
        break;
      default:
        return base::InvalidArgumentError(
            "EXR invalid ValueType. Only Half and Float are supported");
    }
    header.channels().insert(std::string(channel.name),
                             Imf::Channel(pixel_type));

    char* data = nullptr;
    size_t stride = 0;
    if (channel.out_type == ValueType::kHalf) {
      // We want to write half-precision floats to the file, but we are provided
      // with single-precision input.  Perform the conversion here to
      // |half_channel|.
      std::vector<half> half_channel;
      half_channel.reserve(sample_count);
      for (int i = 0; i < sample_count; ++i) {
        half_channel.push_back(channel.data[i * channel.element_period]);
      }
      // Construct |half_channel_ptr| to index into the half-precision data.
      std::vector<const half*> half_channel_ptr;
      half_channel_ptr.reserve(pixel_count);
      int sample_index = 0;
      for (auto sample_count : sample_counts_) {
        half_channel_ptr.push_back(
            sample_count == 0 ? nullptr : half_channel.data() + sample_index);
        sample_index += sample_count;
      }
      data = reinterpret_cast<char*>(half_channel_ptr.data());
      stride = sizeof(half);
      half_channel_list.emplace_back(std::move(half_channel));
      half_channel_ptr_list.emplace_back(std::move(half_channel_ptr));
    } else {
      // We want to write single-precision floats, so just construct
      // |float_channel_ptr| to index into the original data.
      std::vector<const float*> float_channel_ptr;
      float_channel_ptr.reserve(pixel_count);
      int sample_index = 0;
      for (auto sample_count : sample_counts_) {
        int offset = sample_index * channel.element_period;
        float_channel_ptr.push_back(sample_count == 0 ? nullptr
                                                      : channel.data + offset);
        sample_index += sample_count;
      }
      data = reinterpret_cast<char*>(float_channel_ptr.data());
      stride = channel.element_period * sizeof(float);
      float_channel_ptr_list.emplace_back(std::move(float_channel_ptr));
    }

    // Insert the slice for the channel.
    const int x_sampling = 1;
    const int y_sampling = 1;
    const float kDefaultFill = 0.0f;
    framebuffer.insert(std::string(channel.name),
                       Imf::DeepSlice(pixel_type, data, sizeof(void*),
                                      sizeof(void*) * image_width_, stride,
                                      x_sampling, y_sampling, kDefaultFill));
  }

  // Create the EXR file and write out the scanlines.
  image::ExrOStream exr_ostream(exr_contents);
  Imf::DeepScanLineOutputFile file(exr_ostream, header);
  try {
    file.setFrameBuffer(framebuffer);
    file.writePixels(image_height_);
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

DeepExrImageWriter::DeepExrImageWriter() : image_width_(0), image_height_(0) {}

base::Status DeepExrImageWriter::Init(int image_width, int image_height,
                                      std::vector<int> sample_counts) {
  if (image_width <= 0 || image_height <= 0) {
    return base::InvalidArgumentError(absl::Substitute(
        "EXR size ($0, $1) is not positive", image_width, image_height));
  }
  if (sample_counts.size() != image_width * image_height) {
    return base::InvalidArgumentError(absl::Substitute(
        "EXR sample_counts.size()=$0 does not cover image of size ($1, $2)",
        sample_counts.size(), image_width, image_height));
  }
  for (size_t i = 0; i < sample_counts.size(); ++i) {
    if (sample_counts[i] < 0) {
      return base::InvalidArgumentError(absl::Substitute(
          "EXR sample_counts[$0]=$1 is negative", i, sample_counts[i]));
    }
  }
  image_width_ = image_width;
  image_height_ = image_height;
  sample_counts_ = std::move(sample_counts);
  return base::OkStatus();
}

}  // namespace image
}  // namespace seurat
