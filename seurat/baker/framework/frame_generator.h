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

#ifndef VR_SEURAT_BAKER_FRAMEWORK_FRAME_GENERATOR_H_
#define VR_SEURAT_BAKER_FRAMEWORK_FRAME_GENERATOR_H_

#include <algorithm>
#include <array>
#include <string>
#include <vector>

#include "seurat/baker/framework/frame.h"
#include "seurat/baker/framework/frame_sorter.h"
#include "seurat/base/file_system.h"
#include "seurat/ingest/point_cloud_assembler.h"
#include "seurat/tiler/tiler.h"

namespace seurat {
namespace baker {

// Generates Frames for a predefined input.
class FrameGenerator {
 public:
  virtual ~FrameGenerator() = default;

  virtual base::Status GenerateFrames(std::vector<Frame>* frames) = 0;
};

// Wraps a PointCloudAssembler->Tiler->frame pipeline as a FrameGenerator.
class TilerFrameGenerator : public FrameGenerator {
 public:
  // The texture parameterization applied to the tiles.
  enum class TextureParameterization { kObjectSpace = 0, kProjective };

  TilerFrameGenerator(std::unique_ptr<ingest::PointCloudAssembler> merger,
                      std::unique_ptr<tiler::Tiler> tiler,
                      std::unique_ptr<FrameSorter> sorter,
                      TextureParameterization texture_parameterization)
      : merger_(std::move(merger)),
        tiler_(std::move(tiler)),
        sorter_(std::move(sorter)),
        texture_parameterization_(texture_parameterization) {}
  ~TilerFrameGenerator() override = default;

  base::Status GenerateFrames(std::vector<Frame>* frames) override;

 private:
  // Merges the original RGBD images into a simplified point cloud.
  const std::unique_ptr<ingest::PointCloudAssembler> merger_;

  // Generates tiles to approximate the point-sampled geometry.
  const std::unique_ptr<tiler::Tiler> tiler_;

  // Sorts the resulting frames to determine draw order.
  const std::unique_ptr<FrameSorter> sorter_;

  // The texture parameterization applied to the tiles.
  const TextureParameterization texture_parameterization_;
};

// Wraps a FrameGenerator with disk-caching to avoid regenerating frames
// whenever possible.
//
// Tries to load serialized frames using the given FileSystem, or if that fails,
// then runs the fallback and then saves the results to the FileSystem.
//
// Note:  Cache invalidation must be handled externally! i.e. by manually
// deleting the file.
//
// In other words, this is for development/debugging purposes only!
class DiskCachingFrameGenerator : public FrameGenerator {
 public:
  DiskCachingFrameGenerator(std::string filename,
                            std::shared_ptr<base::FileSystem> file_system,
                            std::unique_ptr<FrameGenerator> fallback)
      : filename_(std::move(filename)),
        file_system_(std::move(file_system)),
        fallback_(std::move(fallback)) {}
  ~DiskCachingFrameGenerator() override = default;

  base::Status GenerateFrames(std::vector<Frame>* frames) override;

 private:
  const std::string filename_;
  const std::shared_ptr<base::FileSystem> file_system_;
  const std::unique_ptr<FrameGenerator> fallback_;
};

}  // namespace baker
}  // namespace seurat

#endif  // VR_SEURAT_BAKER_FRAMEWORK_FRAME_GENERATOR_H_
