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

#include "seurat/ingest/pattern_loader.h"

#include <string>

#include "ion/math/matrix.h"
#include "ion/math/matrixutils.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "seurat/base/color.h"
#include "seurat/base/util.h"
#include "seurat/geometry/cube_face.h"

namespace seurat {
namespace ingest {

using base::Color4f;
using geometry::CubeFace;
using ion::math::Matrix2f;
using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Vector2i;

PatternLoader::PatternLoader(const Parameters& parameters)
    : parameters_(parameters) {}

base::Status PatternLoader::LoadViewGroup(
    int view_group_index, std::vector<std::shared_ptr<base::Camera>>* cameras,
    std::vector<image::Ldi4f>* ldis) const {
  const Color4f kWhite(1.0f, 1.0f, 1.0f, 1.0f);
  const Color4f kBlack(0.0f, 0.0f, 0.0f, 1.0f);
  const std::vector<Color4f> kSquareColors = {kWhite, kBlack};
  const Vector2i image_size(parameters_.image_size, parameters_.image_size);
  if (cameras) cameras->clear();
  if (ldis) ldis->clear();
  if (cameras) {
    const Matrix4f clip_from_eye = ion::math::PerspectiveMatrixFromFrustum(
        -1.0f, 1.0f, -1.0f, 1.0f, 1.0f, 100.0f);
    const Matrix4f eye_from_world =
        geometry::LookAtMatrixFromFace(CubeFace::kFront);
    cameras->push_back(std::make_shared<base::ProjectiveCamera>(
        image_size, clip_from_eye, eye_from_world));
  }
  if (ldis) {
    const float angle_radians = parameters_.angle_degrees * M_PI / 180.0f;
    const int feature_size = parameters_.feature_size;
    const float s = std::sin(angle_radians);
    const float c = std::cos(angle_radians);
    const std::string& pattern_name = parameters_.name;
    const Matrix2f plane_rotation(c, -s, s, c);
    std::function<int(const Point2f&)> color_index_function;
    if ("checkerboard" == pattern_name) {
      color_index_function = [=](const Point2f& p) {
        const Point2f rotated_point = plane_rotation * p;
        return static_cast<int>(std::round(rotated_point[0] / feature_size) +
                                std::round(rotated_point[1] / feature_size));
      };
    } else if ("stripes" == pattern_name) {
      color_index_function = [=](const Point2f& p) {
        return static_cast<int>(
            std::round((c * p[0] + s * p[1]) / feature_size));
      };
    } else {
      LOG(FATAL);
    }
    const int num_colors = kSquareColors.size();
    std::vector<Color4f> pixel_colors;
    for (int y = 0; y < image_size[1]; ++y) {
      for (int x = 0; x < image_size[0]; ++x) {
        int color_index = color_index_function(Point2f(x, y));
        if (color_index < 0) {
          // Update |color_index| to a non-negative integer in the same residue
          // class modulo |num_colors|.
          color_index += ((-color_index) / num_colors + 1) * num_colors;
        }
        CHECK_LE(0, num_colors);
        color_index = color_index % num_colors;
        pixel_colors.push_back(kSquareColors[color_index]);
      }
    }
    const int num_pixels = image_size[0] * image_size[1];
    std::vector<int> sample_counts(num_pixels, 1);
    std::vector<float> pixel_depths(num_pixels, parameters_.depth);
    ldis->emplace_back(image_size, std::move(sample_counts),
                       std::move(pixel_colors), std::move(pixel_depths));
  }
  return base::OkStatus();
}

}  // namespace ingest
}  // namespace seurat
