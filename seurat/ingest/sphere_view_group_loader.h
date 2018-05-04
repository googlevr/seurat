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

#ifndef VR_SEURAT_INGEST_SPHERE_VIEW_GROUP_LOADER_H_
#define VR_SEURAT_INGEST_SPHERE_VIEW_GROUP_LOADER_H_

#include <memory>

#include "seurat/base/camera.h"
#include "seurat/image/image.h"
#include "seurat/ingest/view_group_loader.h"

namespace seurat {
namespace ingest {

// Wraps another ViewGroupLoader to overlay rendered spheres.
//
// This only works for single-layer depth images.
class SphereViewGroupLoader : public ViewGroupLoader {
 public:
  struct Sphere {
    ion::math::Point3f center;
    float radius;
  };

  // Determines how to shade a point on a sphere.
  //
  // Example:
  //   SphereShader shader = ...;
  //   Sphere sphere = ...;
  //   Point3f point_on_sphere = ...;
  //   Color3f color = shader(sphere, point_on_sphere);
  using SphereShader =
      std::function<base::Color3f(const Sphere&, const ion::math::Point3f&)>;

  SphereViewGroupLoader(int thread_count,
                        std::unique_ptr<ViewGroupLoader> original,
                        std::vector<Sphere> spheres, SphereShader shader)
      : thread_count_(thread_count),
        original_(std::move(original)),
        spheres_(std::move(spheres)),
        shader_(std::move(shader)) {}

  // ViewGroupLoader implementation.
  int GetNumViewGroups() const override {
    return original_->GetNumViewGroups();
  }
  base::Status LoadViewGroup(
      int view_group_index, std::vector<std::shared_ptr<base::Camera>>* cameras,
      std::vector<image::Ldi4f>* ldis) const override;

 private:
  // Given a camera renders out the spheres to an RGBD image.
  //
  // The returned depth values are stored as Camera::RayDirection scale factors.
  image::Image4f RenderRGBDForCamera(const base::Camera& camera) const;

  // The maximum number of threads to use.
  const int thread_count_;

  // Loads the original ViewGroups to be composited with the rendered spheres.
  const std::unique_ptr<ViewGroupLoader> original_;

  // The spheres to render.
  const std::vector<Sphere> spheres_;

  // Determines how to shade the spheres.
  const SphereShader shader_;
};

}  // namespace ingest
}  // namespace seurat

#endif  // VR_SEURAT_INGEST_SPHERE_VIEW_GROUP_LOADER_H_
