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

#ifndef VR_SEURAT_INGEST_POINT_CLOUD_ASSEMBLER_FACTORY_H_
#define VR_SEURAT_INGEST_POINT_CLOUD_ASSEMBLER_FACTORY_H_

#include <memory>

#include "seurat/ingest/point_cloud_assembler.h"
#include "seurat/ingest/view_group_loader.h"

namespace seurat {
namespace ingest {

class PointCloudAssemblerFactory {
 public:
  struct Parameters {
    // Number of loading threads.
    int thread_count = 1;

    // Resolution of the target display in pixels per degree. This parameter is
    // used to determine binning grid resolution.
    float pixels_per_degree;
  };

  explicit PointCloudAssemblerFactory(const Parameters& parameters)
      : parameters_(parameters) {}

  // Creates the PointCloudAssembler with the given |parameters| and returns it.
  std::unique_ptr<PointCloudAssembler> CreatePointCloudAssembler(
      std::unique_ptr<ViewGroupLoader> view_group_loader);

 private:
  const Parameters parameters_;
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_POINT_CLOUD_ASSEMBLER_FACTORY_H_
