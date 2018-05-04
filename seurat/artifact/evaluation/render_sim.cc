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

#include "seurat/artifact/evaluation/render_sim.h"

#include "ion/math/vector.h"
#include "seurat/base/color.h"
#include "seurat/base/parallel.h"
#include "seurat/geometry/bilinear_interpolator.h"
#include "seurat/geometry/triangle.h"

namespace seurat {
namespace artifact {

using base::Color4f;
using geometry::BilinearInterpolator;
using geometry::Mesh;
using geometry::Raytracer;
using image::Image4f;
using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2f;
using ion::math::Vector2i;
using ion::math::Vector3f;

void RenderSim::Build(const Mesh& mesh, Image4f const* texture) {
  CHECK_EQ(mesh.GetTextureCount(), 1);
  mesh_ = mesh;
  texture_ = texture;
  raytracer_ = Raytracer::Build(mesh_);
}

Color4f RenderSim::TraceRay(const Point3f& origin,
                            const Vector3f& direction) const {
  CHECK(raytracer_) << "RenderSim must be built before it can be used";

  std::vector<Raytracer::Intersection> intersections;
  raytracer_->FindAllIntersections(origin, direction, &intersections);

  if (enable_z_buffer_) {
    // Sort back to front.
    std::sort(intersections.begin(), intersections.end(),
              [](const Raytracer::Intersection& lhs,
                 const Raytracer::Intersection& rhs) {
                return lhs.t_hit > rhs.t_hit;
              });
  } else {
    // Sort by primitive-order.
    //
    // If the primitives are already sorted for rendering with alpha blending
    // via OpenGL, then the result should be ordered back-to-front.
    std::sort(intersections.begin(), intersections.end(),
              [](const Raytracer::Intersection& lhs,
                 const Raytracer::Intersection& rhs) {
                return lhs.triangle_index < rhs.triangle_index;
              });
  }

  // Composite back-to-front using pre-multiplied alpha blending.
  Color4f color = Color4f::Zero();
  for (const auto& intersection : intersections) {
    Point3f point = origin + direction * intersection.t_hit;
    Mesh::Triangle indices = mesh_.GetTriangles()[intersection.triangle_index];
    geometry::Triangle3f triangle_world = {{mesh_.GetPositions()[indices[0]],
                                            mesh_.GetPositions()[indices[1]],
                                            mesh_.GetPositions()[indices[2]]}};

    Point3f bary = geometry::BarycentricFromPoint(triangle_world, point);

    geometry::Triangle3f triangle_texture = {
        {mesh_.GetTexCoords(0)[indices[0]],  //
         mesh_.GetTexCoords(0)[indices[1]],  //
         mesh_.GetTexCoords(0)[indices[2]]}};
    Point3f tex_coords = geometry::PointFromBarycentric(triangle_texture, bary);
    if (!std::isfinite(tex_coords[0]) || !std::isfinite(tex_coords[1]) ||
        !std::isfinite(tex_coords[2])) {
      // If the interpolated texture coordinate ill-defined (e.g. due to grazing
      // angles), then just skip the sample.
      continue;
    }
    Color4f src = SampleTexture(
        Point2f(tex_coords[0] / tex_coords[2], tex_coords[1] / tex_coords[2]));
    Color4f dst = color;
    // Composite 'src over dst' assuming premultiplied alpha textures.
    Color4f out;
    for (int c = 0; c < 4; ++c) {
      out[c] = src[c] + dst[c] * (1.0f - src[3]);
    }
    color = out;
  }

  return color;
}

void RenderSim::Render(const base::Camera& camera,
                       image::Image4f* image) const {
  image->Resize(camera.GetImageSize());
  base::ParallelFor(thread_count_, image->GetSize()[1], [&](int y) {
    for (int x = 0; x < image->GetSize()[0]; ++x) {
      const Point3f origin = camera.RayOrigin({x, y});
      const Vector3f direction = camera.RayDirection({x, y});
      image->At(x, y) = TraceRay(origin, direction);
    }
  });
}

namespace {

void ClampToSize(const Vector2i& size, Point2i* point) {
  (*point)[0] = std::max(0, (*point)[0]);
  (*point)[0] = std::min(size[0] - 1, (*point)[0]);
  (*point)[1] = std::max(0, (*point)[1]);
  (*point)[1] = std::min(size[1] - 1, (*point)[1]);
}

}  // namespace

Color4f RenderSim::SampleTexture(const Point2f& texture_coordinate) const {
  Vector2i size = texture_->GetSize();
  // Pixel coordinate relative to the texture's integer coordinates.
  //
  // In other words, p = (0,0) corresponds exactly to the sample retrieved via
  // 'texture_->At(0,0)'.
  //
  // This can be seen as converting from OpenGL texture coordinates to a
  // coordinate space relative to the pixels of the image.
  Point2f p(texture_coordinate[0] * size[0] - 0.5f,
            texture_coordinate[1] * size[1] - 0.5f);

  // Integer coordinates of the 4 nearest pixels to use for linear resampling.
  Point2i s00(p[0], p[1]);
  Point2i s10(p[0] + 1.0f, p[1]);
  Point2i s01(p[0], p[1] + 1.0f);
  Point2i s11(p[0] + 1.0f, p[1] + 1.0f);

  // Clamp to the edge (i.e. GL_CLAMP_TO_EDGE instead of GL_REPEAT).
  ClampToSize(size, &s00);
  ClampToSize(size, &s10);
  ClampToSize(size, &s01);
  ClampToSize(size, &s11);

  BilinearInterpolator<Color4f> interp(
      std::array<Color4f, 4>{{texture_->At(s00), texture_->At(s10),
                              texture_->At(s11), texture_->At(s01)}});

  return interp.At(p - Point2f(Point2i(p)) + Point2f::Zero());
}

}  // namespace artifact
}  // namespace seurat
