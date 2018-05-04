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

#ifndef VR_SEURAT_INGEST_VIEW_GROUP_LOADER_TEST_UTILS_H_
#define VR_SEURAT_INGEST_VIEW_GROUP_LOADER_TEST_UTILS_H_

#include <array>
#include <memory>

#include "ion/math/vector.h"
#include "seurat/base/color.h"
#include "seurat/ingest/view_group_loader.h"

namespace seurat {
namespace ingest {

// A view group loader for testing. It returns six faces of a cube map for
// each view group. The color and depth images are set to uniform values. The
// values for each face can be configured.
class FakeViewGroupLoader : public ViewGroupLoader {
 public:
  // Creates a view group loader that returns six faces of a cube map for each
  // of the |num_view_groups| view groups. |image_size| specifies the resolution
  // of each cube face. The color and depth values of each face are defined by
  // |face_colors| and |face_depths|, respectively.
  FakeViewGroupLoader(int num_view_groups,
                      const ion::math::Vector2i& image_size,
                      const std::array<base::Color3f, 6>& face_colors,
                      const std::array<float, 6>& face_depths);

  // Creates a view group loader for testing cameras. It returns six faces of a
  // cube map for each view groups. The world-space positions of the view groups
  // are passed in as |positions|. |image_size| specifies the resolution of each
  // cube face. Colors and depths are set to black and 0, respectively.
  FakeViewGroupLoader(const std::vector<ion::math::Point3f>& positions,
                      const ion::math::Vector2i& image_size);

  // ViewGroupLoader implementation.
  int GetNumViewGroups() const override { return num_view_groups_; }
  base::Status LoadViewGroup(
      int view_group_index, std::vector<std::shared_ptr<base::Camera>>* cameras,
      std::vector<image::Ldi4f>* ldis) const override;

 private:
  // The number of view groups returned by this loader.
  const int num_view_groups_;

  // World-space positions of the view groups.
  const std::vector<ion::math::Point3f> positions_;

  // The resolution of each cube face.
  const ion::math::Vector2i image_size_;

  // The color value for each of the six cube faces. The color inside a face is
  // uniform.
  const std::array<base::Color3f, 6> face_colors_;

  // The depth value for each of the six cube faces. The depth inside a face is
  // uniform.
  const std::array<float, 6> face_depths_;
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_VIEW_GROUP_LOADER_TEST_UTILS_H_
