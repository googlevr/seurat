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

#ifndef VR_SEURAT_BAKER_FRAMEWORK_RAY_BUNDLE_H_
#define VR_SEURAT_BAKER_FRAMEWORK_RAY_BUNDLE_H_

#include <tuple>
#include <vector>

#include "ion/math/vector.h"
#include "absl/types/span.h"
#include "seurat/base/camera.h"
#include "seurat/base/color.h"
#include "seurat/image/ldi.h"
#include "seurat/ingest/view_group_loader.h"

namespace seurat {
namespace baker {

// A generalized layered depth image.
//
// The interface presents a structure-of-arrays interface in which each Ray
// (specified by an index) has:
//  * An origin.
//  * A set of intersections with scene geometry (aka "samples"), sorted in
//    increasing distance.
//
// And each intersection has a point & RGBA value.
//
// For example, an implementation may be backed by one of:
//  * A single RGB+D image.
//  * An RGB+D omni-directional stereo panorama.
//  * A layered depth image.
//  * A cubemap of layered depth images.
//  * Incoherent/unstructured rays.
class RayBundle {
 public:
  // A (ray index, intersection index) tuple.
  using RayIntersectionIndex = std::tuple<int, int>;

  virtual ~RayBundle() = default;

  // Returns the total number of rays.
  virtual int GetRayCount() const = 0;

  // Returns the origin of the specified ray.
  virtual ion::math::Point3f GetOrigin(int ray_index) const = 0;

  // Returns the direction of the specified ray.
  virtual ion::math::Vector3f GetDirection(int ray_index) const = 0;

  // Returns the number of intersections of the specified ray with scene
  // geometry.
  //
  // For example, an implementation backed by a single depth image will always
  // return 0-1.  A layered-depth-image implementation may return [0, 1, ..., ].
  virtual int GetIntersectionCount(int ray_index) const = 0;

  // Returns the point of intersection of the specified ray.
  //
  // All endpoints for a given ray *must* be colinear.
  //
  // |intersection_index| must be in [0, GetIntersectionCount(ray_index) - 1].
  virtual ion::math::Point3f GetIntersectionPoint(
      int ray_index, int intersection_index) const = 0;

  // Returns the observed color at the specified intersection of the ray.
  //
  // |intersection_index| must be in [0, GetIntersectionCount(ray_index) - 1].
  virtual base::Color4f GetIntersectionColor(int ray_index,
                                             int intersection_index) const = 0;
};

// Wraps a ViewGroupLoader to load into a RayBundle.
class ViewGroupRayBundle : public RayBundle {
 public:
  ViewGroupRayBundle() = default;
  ViewGroupRayBundle(std::vector<std::shared_ptr<base::Camera>> cameras,
                     std::vector<image::Ldi4f> ldis);

  ~ViewGroupRayBundle() override = default;

  // RayBundle implementation.
  int GetRayCount() const override;
  ion::math::Point3f GetOrigin(int ray_index) const override;
  ion::math::Vector3f GetDirection(int ray_index) const override;
  int GetIntersectionCount(int ray_index) const override;
  ion::math::Point3f GetIntersectionPoint(
      int ray_index, int intersection_index) const override;
  base::Color4f GetIntersectionColor(int ray_index,
                                     int intersection_index) const override;

 private:
  // Returns the index of the view & pixel-coordinates corresponding to the
  // specified ray.
  std::tuple<int, ion::math::Point2i> PixelCoordFromRayIndex(
      int ray_index) const;

  // The ViewGroup:  A set of parallel Cameras and Ldis.
  std::vector<std::shared_ptr<base::Camera>> cameras_;
  std::vector<image::Ldi4f> ldis_;

  ion::math::Vector2i image_size_;
  int num_views_;
};

// The simplest RayBundle which explicitly stores a collection of rays
// internally in a structure-of-arrays format.
class ExplicitRayBundle : public RayBundle {
 public:
  ExplicitRayBundle() = default;
  ~ExplicitRayBundle() override = default;

  // Appends a ray.
  void AddRay(const ion::math::Point3f& origin,
              const ion::math::Vector3f& direction,
              absl::Span<const float> intersection_depths,
              absl::Span<const base::Color4f> intersection_colors) {
    CHECK_EQ(intersection_depths.size(), intersection_colors.size());
    origins_.push_back(origin);
    directions_.push_back(direction);
    intersection_depths_.emplace_back(intersection_depths.begin(),
                                      intersection_depths.end());
    intersection_colors_.emplace_back(intersection_colors.begin(),
                                      intersection_colors.end());
  }

  // RayBundle implementation.
  int GetRayCount() const override { return origins_.size(); }
  ion::math::Point3f GetOrigin(int ray_index) const override {
    return origins_[ray_index];
  }
  ion::math::Vector3f GetDirection(int ray_index) const override {
    return directions_[ray_index];
  }
  int GetIntersectionCount(int ray_index) const override {
    return intersection_depths_[ray_index].size();
  }
  ion::math::Point3f GetIntersectionPoint(
      int ray_index, int intersection_index) const override {
    return intersection_depths_[ray_index][intersection_index] *
               directions_[ray_index] +
           origins_[ray_index];
  }
  base::Color4f GetIntersectionColor(int ray_index,
                                     int intersection_index) const override {
    return intersection_colors_[ray_index][intersection_index];
  }

 private:
  // The origin of each ray.
  std::vector<ion::math::Point3f> origins_;
  // The direction of each ray.
  std::vector<ion::math::Vector3f> directions_;
  // The depths of each intersection (as a multiple of directions_) of each ray.
  std::vector<std::vector<float>> intersection_depths_;
  // The color of the point of each intersection.
  //
  // This runs parallel to intersection_depths_.
  std::vector<std::vector<base::Color4f>> intersection_colors_;
};

}  // namespace baker
}  // namespace seurat

#endif  // VR_SEURAT_BAKER_FRAMEWORK_RAY_BUNDLE_H_
