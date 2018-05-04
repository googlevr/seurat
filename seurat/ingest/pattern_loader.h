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

#ifndef VR_SEURAT_INGEST_PATTERN_LOADER_H_
#define VR_SEURAT_INGEST_PATTERN_LOADER_H_

#include <string>

#include "ion/math/vector.h"
#include "seurat/ingest/view_group_loader.h"

namespace seurat {
namespace ingest {

// A view group loader for testing. It returns a single view group representing
// the front face of a cube map. The colors are set to a checkerboard pattern
// and the depth is set to a constant.
class PatternLoader : public ViewGroupLoader {
 public:
  // This struct wraps the parameters that control the pattern.
  struct Parameters {
    // Pattern name. Currently one of "checkerboard" or "stripes".
    std::string name;

    // LDI resolution.
    int image_size;

    // The size of a checkerboard square, or of a stripe, in units of image
    // pixels.
    int feature_size;

    // Pattern's angle relative to the LDI's raster lines, in degrees.
    float angle_degrees;

    // Z-depth of LDI samples.
    float depth;
  };

  explicit PatternLoader(const Parameters& parameters);

  int GetNumViewGroups() const override { return 1; }

  base::Status LoadViewGroup(
      int view_group_index, std::vector<std::shared_ptr<base::Camera>>* cameras,
      std::vector<image::Ldi4f>* ldis) const override;

 private:
  // Parameters controlling the pattern.
  const Parameters parameters_;
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_PATTERN_LOADER_H_
