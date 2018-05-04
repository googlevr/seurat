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

#ifndef VR_SEURAT_INGEST_VIEW_GROUP_LOADER_FACTORY_H_
#define VR_SEURAT_INGEST_VIEW_GROUP_LOADER_FACTORY_H_

#include <limits>
#include <memory>
#include <string>

#include "absl/strings/string_view.h"
#include "seurat/base/file_system.h"
#include "seurat/base/status.h"
#include "seurat/ingest/pattern_loader.h"
#include "seurat/ingest/view_group_loader.h"

namespace seurat {
namespace ingest {

class ViewGroupLoaderFactory {
 public:
  struct Parameters {
    // Number of loading threads.
    int thread_count = 1;

    // If not empty, process only the specified face. Must be one of 'front',
    // 'back', 'left', 'right', 'bottom', 'top'.
    std::string single_face = "";

    // Maximum number of views to load.
    int max_views = std::numeric_limits<int>::max();

    // Half the side-length of the origin-centered skybox to clamp geometry.
    //
    // 0 indicates no skybox clamping should be performed.
    float skybox_radius = 0.0f;

    // Parameters for an analytic pattern, if any.
    PatternLoader::Parameters pattern_parameters;
  };

  ViewGroupLoaderFactory(const Parameters& parameters,
                         std::shared_ptr<base::FileSystem> file_system)
      : parameters_(parameters), file_system_(std::move(file_system)) {}

  // Creates the view group loader with the given json or proto
  // |manifest_filename| and returns it in |view_group_loader|.
  base::Status CreateViewGroupLoader(
      absl::string_view manifest_filename,
      std::unique_ptr<ViewGroupLoader>* view_group_loader) const;

 private:
  // Parameters applied to all constructed loaders.
  const Parameters parameters_;

  // The filesystem used to load manifests & other assets.
  const std::shared_ptr<base::FileSystem> file_system_;
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_VIEW_GROUP_LOADER_FACTORY_H_
