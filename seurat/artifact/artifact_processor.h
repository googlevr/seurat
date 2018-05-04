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

#ifndef VR_SEURAT_ARTIFACT_ARTIFACT_PROCESSOR_H_
#define VR_SEURAT_ARTIFACT_ARTIFACT_PROCESSOR_H_

#include <memory>
#include <vector>

#include "seurat/artifact/artifact.h"
#include "seurat/base/status.h"
#include "seurat/base/status_util.h"

namespace seurat {
namespace artifact {

class ArtifactProcessor {
 public:
  virtual ~ArtifactProcessor() = default;

  // Processes an artifact.
  virtual base::Status Process(Artifact* artifact) const {
    return base::OkStatus();
  }
};

// In the ASCII art documentation below, "-" denotes a sequence and "/" and "\"
// denote a group. "{foo}" denotes an artifact. "node_xyz" denotes an artifact
// processor. Arrows "->" denote input and output data flow.

// A sequence of artifact processors (sequential pipeline). The output of one
// stage is fed as input into the next stage. ArtifactProcessorSequence is an
// artifact transform.
//
// Sequence data flow:
//
// {input}->node_0->{inout_0}->node_1->{inout_1}->...->node_N-1->{output}
//
class ArtifactProcessorSequence : public ArtifactProcessor {
 public:
  explicit ArtifactProcessorSequence(
      std::vector<std::shared_ptr<const ArtifactProcessor>> stages)
      : stages_(std::move(stages)) {}

  base::Status Process(Artifact* artifact) const override {
    for (const auto& stage : stages_) {
      SEURAT_RETURN_IF_ERROR(stage->Process(artifact));
    }
    return base::OkStatus();
  }

 private:
  // The stages of the sequence.
  std::vector<std::shared_ptr<const ArtifactProcessor>> stages_;
};

// A group of artifact processors (branching pipeline). ArtifactProcessorGroup
// is an artifact sink.
//
// Group data flow:
//
//         >node_0
//        /
// {input}
//        \
//         >node_1
//
// Pipeline example that combines groups and sequences:
//
//                    node_0
//                   /
//     node_3-node_2-
//    /              \
//   / node_4         node_1
//   \/
//    \
//     node_5-node_6
//
class ArtifactProcessorGroup : public ArtifactProcessor {
 public:
  explicit ArtifactProcessorGroup(
      std::vector<std::shared_ptr<const ArtifactProcessor>> children)
      : children_(std::move(children)) {}

  base::Status Process(Artifact* artifact) const override {
    base::Status group_status;
    for (const auto& child : children_) {
      base::Status child_status;
      Artifact copy = *artifact;
      child_status = child->Process(&copy);
      base::UpdateStatus(&group_status, child_status);
    }
    return group_status;
  }

 private:
  // That children of this node.
  std::vector<std::shared_ptr<const ArtifactProcessor>> children_;
};

}  // namespace artifact
}  // namespace seurat

#endif  // VR_SEURAT_ARTIFACT_ARTIFACT_PROCESSOR_H_
