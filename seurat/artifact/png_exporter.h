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

#ifndef VR_SEURAT_ARTIFACT_PNG_EXPORTER_H_
#define VR_SEURAT_ARTIFACT_PNG_EXPORTER_H_

#include <memory>
#include <string>

#include "seurat/artifact/artifact.h"
#include "seurat/artifact/artifact_processor.h"
#include "seurat/base/file_system.h"

namespace seurat {
namespace artifact {

// Writes Artifacts to disk as .png files storing their texture.
class PngExporter : public ArtifactProcessor {
 public:
  explicit PngExporter(std::shared_ptr<base::FileSystem> filesystem,
                       std::string basename)
      : filesystem_(std::move(filesystem)), basename_(std::move(basename)) {}
  ~PngExporter() override = default;

  // Exporter implementation.
  base::Status Process(Artifact* artifact) const override;

 private:
  // Used to write to disk.
  const std::shared_ptr<base::FileSystem> filesystem_;

  // The path + basename to the output image.  A file extension will be added by
  // the artifact to determine the complete path/name.
  const std::string basename_;
};

}  // namespace artifact
}  // namespace seurat

#endif  // VR_SEURAT_ARTIFACT_PNG_EXPORTER_H_
