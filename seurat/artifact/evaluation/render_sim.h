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

#ifndef VR_SEURAT_ARTIFACT_EVALUATION_RENDER_SIM_H_
#define VR_SEURAT_ARTIFACT_EVALUATION_RENDER_SIM_H_

#include <memory>

#include "ion/math/vector.h"
#include "seurat/base/camera.h"
#include "seurat/base/color.h"
#include "seurat/geometry/raytracer.h"
#include "seurat/image/image.h"

namespace seurat {
namespace artifact {

// OpenGL Rendering Simulator for typical Seurat output (textured meshes).
//
// A raytracer which simulates rendering with OpenGL, complete with
// alpha-sorting artifacts (i.e. rendering based on *primitive order* instead of
// depth) and poor-quality bilinear texture filtering.
//
// This is intended for use in end-to-end benchmarking, evaluation, and testing.
class RenderSim {
 public:
  explicit RenderSim(const int thread_count, bool enable_z_buffer = false)
      : enable_z_buffer_(enable_z_buffer), thread_count_(thread_count) {}

  // Initializes with the given |mesh| and |texture|, building any acceleration
  // structures required for rendering.
  //
  // |mesh| must have texture coordinates pointing into the |texture| using the
  // same conventions as OpenGL.
  //
  // |texture| is assumed to be premultiplied alpha, as would be expected if
  // rendering using OpenGL to get ideal texture filtering.
  //
  // |texture| must persist through all subsequent calls to TraceRay().
  void Build(const geometry::Mesh& mesh, image::Image4f const* texture);

  // Renders the specified ray.
  //
  // The renderer must be built before this can be used.
  base::Color4f TraceRay(const ion::math::Point3f& ray_origin,
                         const ion::math::Vector3f& ray_direction) const;

  // Renders an image for the given |camera| view.
  void Render(const base::Camera& camera, image::Image4f* image) const;

 private:
  // Samples from texture_, modeling OpenGL's nearest-neighbor resampling.
  base::Color4f SampleTexture(
      const ion::math::Point2f& texture_coordinate) const;

  // If true, renders as if there is a z-buffer.
  //
  // If false, renders based on the order of the primitives.  This simulates the
  // expected behavior when rendering alpha-blended geometry with OpenGL on a
  // mobile device.
  const bool enable_z_buffer_;

  // The number of threads to use when rendering.
  const int thread_count_;

  // The mesh currently being rendered.
  geometry::Mesh mesh_;

  // The texture associated with the mesh currently being rendered.
  image::Image4f const* texture_;

  // Traces rays through the triangles in mesh_.
  std::unique_ptr<geometry::Raytracer> raytracer_;
};

}  // namespace artifact
}  // namespace seurat

#endif  // VR_SEURAT_ARTIFACT_EVALUATION_RENDER_SIM_H_
