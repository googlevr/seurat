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

#include "seurat/ingest/single_face_view_group_loader.h"

#include "ion/math/matrix.h"
#include "ion/math/matrixutils.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "seurat/base/camera.h"
#include "seurat/base/color.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/base/parallel.h"
#include "seurat/base/status.h"

namespace seurat {
namespace ingest {

using base::Color4f;
using image::Ldi4f;
using ion::math::Point3f;
using ion::math::Vector3f;

base::Status SingleFaceViewGroupLoader::LoadViewGroup(
    int view_group_index, std::vector<std::shared_ptr<base::Camera>>* cameras,
    std::vector<image::Ldi4f>* ldis) const {
  std::vector<std::shared_ptr<base::Camera>> camera_buffer;

  base::Status status =
      delegate_->LoadViewGroup(view_group_index, &camera_buffer, ldis);

  const Point3f cube_face_forward =
      ion::math::Inverse(geometry::LookAtMatrixFromFace(cube_face_)) *
      Point3f(0.0f, 0.0f, -1.0f);

  const int major_axis = base::MajorAxisFromPosition(cube_face_forward);
  const int axis1 = (major_axis + 1) % 3;
  const int axis2 = (major_axis + 2) % 3;
  const int sign = cube_face_forward[major_axis] > 0 ? 1 : -1;

  if (ldis) {
    base::ParallelFor(thread_count_, ldis->size(), [&](int i) {
      const auto& camera = camera_buffer[i];
      const Ldi4f& ldi = ldis->at(i);
      const ion::math::Vector2i size = ldi.GetSize();
      std::vector<int> sample_counts;
      std::vector<Color4f> colors;
      std::vector<float> depths;
      for (int y = 0; y < ldi.GetHeight(); ++y) {
        for (int x = 0; x < ldi.GetWidth(); ++x) {
          int samples_added = 0;
          for (int s = 0; s < ldi.GetSampleCount({x, y}); ++s) {
            Point3f p = camera->RayEnd({x, y}, ldi.GetDepths({x, y})[s]);
            float depth = p[major_axis] * sign;
            if (depth > 0 && depth > std::fabs(p[axis1]) &&
                depth > std::fabs(p[axis2])) {
              colors.push_back(ldi.GetColors({x, y})[s]);
              depths.push_back(ldi.GetDepths({x, y})[s]);
              samples_added++;
            }
          }
          sample_counts.push_back(samples_added);
        }
      }
      ldis->at(i) = image::Ldi4f(size, sample_counts, std::move(colors),
                                 std::move(depths));
    });
  }

  if (cameras) {
    *cameras = std::move(camera_buffer);
  }

  return status;
}

}  // namespace ingest
}  // namespace seurat
