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

#include "seurat/baker/framework/frame_generator.h"

#include "ion/math/vector.h"
#include "seurat/base/file_system.h"
#include "seurat/base/parallel.h"
#include "seurat/base/reporting.h"
#include "seurat/base/util.h"
#include "seurat/geometry/cube_face.h"
#include "seurat/tiler/point_set.h"

namespace seurat {
namespace baker {

using base::FileSystem;
using ion::math::Point3f;
using tiler::PointSet;

base::Status TilerFrameGenerator::GenerateFrames(std::vector<Frame>* frames) {
  std::vector<ion::math::Point3f> positions;
  std::vector<float> weights;
  SEURAT_RETURN_IF_ERROR(merger_->Build(&positions, &weights));

  PointSet point_set;
  point_set.id = 0;
  point_set.positions = positions;
  point_set.weights = weights;
  std::vector<tiler::Tile> tiles = tiler_->Run(point_set);

  for (const auto& tile : tiles) {
    Frame frame;
    frame.draw_order = -1;
    frame.quad = tile.quad;
    switch (texture_parameterization_) {
      case TextureParameterization::kProjective: {
        const Point3f frame_center =
            0.25f * (tile.quad[0] + tile.quad[1] + tile.quad[2] + tile.quad[3]);
        const geometry::CubeFace cube_face =
            geometry::CubeFaceFromPosition(frame_center);
        for (int i = 0; i < 4; ++i) {
          const Point3f uvw = geometry::UVWFromXYZ(tile.quad[i], cube_face);
          frame.texcoord_w[i] = uvw[2];
        }
        break;
      }
      case TextureParameterization::kObjectSpace: {
        std::fill(frame.texcoord_w.begin(), frame.texcoord_w.end(), 1.0f);
        break;
      }
    }
    frames->push_back(frame);
  }

  sorter_->ComputeDrawOrder(positions, absl::MakeSpan(*frames));

  return base::OkStatus();
}

namespace {

// Reads a vector of packed Frames from the given string, assumed to be valid.
std::vector<Frame> FramesFromString(const std::string& cached_frames) {
  std::vector<Frame> frames(cached_frames.size() / sizeof(Frame));
  for (int frame_index = 0; frame_index < frames.size(); ++frame_index) {
    int offset = sizeof(Frame) * frame_index;
    frames[frame_index] =
        *reinterpret_cast<const Frame*>(&cached_frames[offset]);
  }
  return frames;
}

// Serializes a vector of frames to a string.
std::string SerializeFramesToString(const std::vector<Frame>& frames) {
  // The initialization character ('z') is irrelevant, since it should all be
  // overwritten.  Using a non-zero character may help if debugging.
  std::string encoded_frames(frames.size() * sizeof(Frame), 'z');

  for (int frame_index = 0; frame_index < frames.size(); ++frame_index) {
    int offset = sizeof(Frame) * frame_index;
    const char* start = reinterpret_cast<const char*>(&(frames[frame_index]));
    const char* end = start + sizeof(Frame);
    std::copy(start, end, &encoded_frames[offset]);
  }
  return encoded_frames;
}

}  // namespace

base::Status DiskCachingFrameGenerator::GenerateFrames(
    std::vector<Frame>* frames) {
  // Load & return cached frames.
  std::string cache_contents;
  if (file_system_->GetContents(filename_, &cache_contents).ok()) {
    if (cache_contents.size() % sizeof(Frame) != 0) {
      base::SeuratWarning("Invalid geometry cache: " +
                             file_system_->GetAbsolutePath(filename_));
    } else {
      *frames = FramesFromString(cache_contents);
      base::SeuratInfo(
          "Loaded geometry from file: " +
          file_system_->GetAbsolutePath(filename_) +
          ".  Number of quads: " + std::to_string(frames->size()));
      return base::OkStatus();
    }
  }

  // If we get here, then we failed to load cached frames.
  base::SeuratInfo("Unable to load geometry from cache: " +
                      file_system_->GetAbsolutePath(filename_));
  base::SeuratInfo("Recomputing geometry");
  SEURAT_RETURN_IF_ERROR(fallback_->GenerateFrames(frames));
  SEURAT_RETURN_IF_ERROR(
      file_system_->SetContents(filename_, SerializeFramesToString(*frames)));
  base::SeuratInfo("Saved geometry to cache: " +
                      file_system_->GetAbsolutePath(filename_));

  return base::OkStatus();
}

}  // namespace baker
}  // namespace seurat
