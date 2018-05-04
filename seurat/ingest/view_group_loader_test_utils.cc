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

#include "seurat/ingest/view_group_loader_test_utils.h"

#include <vector>

#include "ion/math/matrix.h"
#include "ion/math/transformutils.h"
#include "seurat/base/projective_camera.h"
#include "seurat/base/status.h"
#include "seurat/base/util.h"
#include "seurat/geometry/cube_face.h"
#include "seurat/image/image.h"

namespace seurat {
namespace ingest {

using base::Color3f;
using base::Color4f;
using geometry::CubeFace;
using ion::math::Matrix4f;
using ion::math::Point3f;
using ion::math::Vector2i;

FakeViewGroupLoader::FakeViewGroupLoader(
    int num_view_groups, const Vector2i& image_size,
    const std::array<Color3f, 6>& face_colors,
    const std::array<float, 6>& face_depths)
    : num_view_groups_(num_view_groups),
      positions_(num_view_groups, Point3f::Zero()),
      image_size_(image_size),
      face_colors_(face_colors),
      face_depths_(face_depths) {}

FakeViewGroupLoader::FakeViewGroupLoader(const std::vector<Point3f>& positions,
                                         const Vector2i& image_size)
    : num_view_groups_(positions.size()),
      positions_(positions),
      image_size_(image_size),
      face_colors_({{Color3f::Zero(), Color3f::Zero(), Color3f::Zero(),
                     Color3f::Zero(), Color3f::Zero(), Color3f::Zero()}}),
      face_depths_({{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}}) {}

base::Status FakeViewGroupLoader::LoadViewGroup(
    int view_group_index, std::vector<std::shared_ptr<base::Camera>>* cameras,
    std::vector<image::Ldi4f>* ldis) const {
  // Populate the view group with six faces of a cube map.
  if (cameras) cameras->clear();
  if (ldis) ldis->clear();
  for (int face_index = 0; face_index < 6; ++face_index) {
    const CubeFace cube_face = static_cast<CubeFace>(face_index);
    if (cameras) {
      const Matrix4f clip_from_eye = ion::math::PerspectiveMatrixFromFrustum(
          -1.0f, 1.0f, -1.0f, 1.0f, 1.0f, 100.0f);
      const Matrix4f eye_from_world =
          geometry::LookAtMatrixFromFace(cube_face) *
          ion::math::TranslationMatrix(Point3f::Zero() -
                                       positions_[view_group_index]);
      cameras->push_back(std::make_shared<base::ProjectiveCamera>(
          image_size_, clip_from_eye, eye_from_world));
    }
    // Fill the color image and depth image with the uniform values for the
    // current face.
    if (ldis) {
      const int num_pixels = image_size_[0] * image_size_[1];
      std::vector<int> sample_counts(num_pixels, 1);
      std::vector<Color4f> face_colors(num_pixels,
                                       Color4f(face_colors_[face_index], 1.0f));
      std::vector<float> face_depths(num_pixels, face_depths_[face_index]);
      ldis->emplace_back(image_size_, std::move(sample_counts),
                         std::move(face_colors), std::move(face_depths));
    }
  }
  return base::OkStatus();
}

}  // namespace ingest
}  // namespace seurat
