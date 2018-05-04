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

#include "seurat/baker/framework/ray_bundle.h"

#include "ion/math/transformutils.h"
#include "seurat/base/status.h"
#include "seurat/base/util.h"

namespace seurat {
namespace baker {

using base::Color4f;
using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2f;
using ion::math::Vector2f;
using ion::math::Vector3f;

ViewGroupRayBundle::ViewGroupRayBundle(
    std::vector<std::shared_ptr<base::Camera>> cameras,
    std::vector<image::Ldi4f> ldis)
    : cameras_(std::move(cameras)), ldis_(std::move(ldis)) {
  CHECK_EQ(cameras_.size(), ldis_.size());

  num_views_ = cameras_.size();
  if (num_views_ == 0) {
    image_size_ = {0, 0};
  } else {
    image_size_ = cameras_.front()->GetImageSize();
    for (int i = 0; i < num_views_; ++i) {
      CHECK_EQ(image_size_, cameras_[i]->GetImageSize());
      CHECK_EQ(image_size_, ldis_[i].GetSize());
    }
  }
}

int ViewGroupRayBundle::GetRayCount() const {
  return image_size_[0] * image_size_[1] * num_views_;
}

Point3f ViewGroupRayBundle::GetOrigin(int ray_index) const {
  int view_index;
  Point2i pixel_coords;
  std::tie(view_index, pixel_coords) = PixelCoordFromRayIndex(ray_index);
  return cameras_[view_index]->RayOrigin(pixel_coords);
}

Vector3f ViewGroupRayBundle::GetDirection(int ray_index) const {
  // Rays originate from the near clipping plane.
  int view_index;
  Point2i pixel_coords;
  std::tie(view_index, pixel_coords) = PixelCoordFromRayIndex(ray_index);
  return cameras_[view_index]->RayDirection(pixel_coords);
}

int ViewGroupRayBundle::GetIntersectionCount(int ray_index) const {
  int view_index;
  Point2i pixel_coords;
  std::tie(view_index, pixel_coords) = PixelCoordFromRayIndex(ray_index);
  return ldis_[view_index].GetSampleCount(pixel_coords);
}

Point3f ViewGroupRayBundle::GetIntersectionPoint(int ray_index,
                                                 int intersection) const {
  int view_index;
  Point2i pixel_coords;
  std::tie(view_index, pixel_coords) = PixelCoordFromRayIndex(ray_index);
  float depth = ldis_[view_index].GetDepths(pixel_coords)[intersection];
  return cameras_[view_index]->RayEnd(pixel_coords, depth);
}

Color4f ViewGroupRayBundle::GetIntersectionColor(int ray_index,
                                                 int intersection) const {
  int view_index;
  Point2i pixel_coords;
  std::tie(view_index, pixel_coords) = PixelCoordFromRayIndex(ray_index);
  return ldis_[view_index].GetColors(pixel_coords)[intersection];
}

std::tuple<int, Point2i> ViewGroupRayBundle::PixelCoordFromRayIndex(
    int ray_index) const {
  int rays_per_image = image_size_[0] * image_size_[1];
  int ray_index_in_image = ray_index % rays_per_image;
  return std::make_tuple(ray_index / rays_per_image,
                         Point2i(ray_index_in_image % image_size_[0],
                                 ray_index_in_image / image_size_[0]));
}

}  // namespace baker
}  // namespace seurat
