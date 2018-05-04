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

#ifndef VR_SEURAT_PIPELINE_PIPELINE_H_
#define VR_SEURAT_PIPELINE_PIPELINE_H_

#include <functional>
#include <limits>
#include <memory>
#include <string>

#include "seurat/artifact/artifact.h"
#include "seurat/artifact/artifact_processor.h"
#include "seurat/pipeline/flags.pb.h"

namespace seurat {
namespace pipeline {

// Assembles the processing pipeline.
class Pipeline {
 public:
  // Exposes functions for lazily-constructing various entry points into
  // Seurat.
  struct Factories {
    // Constructs the ProcessingPipeline.
    std::function<std::shared_ptr<const artifact::ArtifactProcessor>()>
        processing_pipeline;
  };

  explicit Pipeline(const proto::Flags& f) : flags_(f) {}

  // Wires together factory functions which lazily-construct pieces of the
  // pipeline.
  Factories AssemblePipelineFactories() const;

 private:
  const proto::Flags flags_;
};

}  // namespace pipeline
}  // namespace seurat

#endif  // VR_SEURAT_PIPELINE_PIPELINE_H_
