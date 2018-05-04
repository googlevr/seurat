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

#include "seurat/ingest/sphere_view_group_loader.h"

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
#include "seurat/geometry/ray_sphere_intersection.h"

namespace seurat {
namespace ingest {

using base::Camera;
using base::Color3f;
using base::Color4f;
using image::Image4f;
using image::Ldi4f;
using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2i;
using ion::math::Vector3f;

namespace {

// Wraps any other Camera type to treat "depth" as "Ray Depth".
//
// Formally, the 'depth' parameter passed in to RayEndFloat should be of the
// form 't' in:
//   RayEnd = RayOrigin + t * RayDirection
class RayDepthCamera : public base::Camera {
 public:
  explicit RayDepthCamera(std::shared_ptr<Camera> delegate)
      : delegate_(std::move(delegate)) {}

  // Camera implementation.
  //
  // RayEndFloat() is the only nontrivial implementation.
  Point3f RayEndFloat(const Point2f& pixel, float depth) const override {
    return RayOriginFloat(pixel) + RayDirectionFloat(pixel) * depth;
  }
  Vector2i GetImageSize() const override { return delegate_->GetImageSize(); }
  Matrix4f GetWorldFromEye() const override {
    return delegate_->GetWorldFromEye();
  }
  Point3f RayOriginFloat(const Point2f& pixel) const override {
    return delegate_->RayOriginFloat(pixel);
  }
  Vector3f RayDirectionFloat(const Point2f& pixel) const override {
    return delegate_->RayDirectionFloat(pixel);
  }

 private:
  const std::shared_ptr<Camera> delegate_;
};

// Converts depth values in the |ldi| to depths compatible with RayDepthCamera.
Ldi4f ConvertToRayDepth(const Ldi4f& ldi, const Camera& original_camera) {
  std::vector<int> sample_counts;
  std::vector<Color4f> colors;
  std::vector<float> depths;
  for (int y = 0; y < ldi.GetHeight(); ++y) {
    for (int x = 0; x < ldi.GetWidth(); ++x) {
      // Pass through the original sample-counts & colors.
      sample_counts.push_back(ldi.GetSampleCount({x, y}));
      for (int s = 0; s < ldi.GetSampleCount({x, y}); ++s) {
        colors.push_back(ldi.GetColors({x, y})[s]);
      }

      // Convert depths to correspond to 't' in the equation:
      //   RayEnd = RayOrigin + t * RayDirection
      const Point3f ray_origin = original_camera.RayOrigin({x, y});
      const Vector3f ray_direction = original_camera.RayDirection({x, y});
      float ray_direction_length2 = ion::math::LengthSquared(ray_direction);
      for (int s = 0; s < ldi.GetSampleCount({x, y}); ++s) {
        Point3f ray_end =
            original_camera.RayEnd({x, y}, ldi.GetDepths({x, y})[s]);
        float new_depth = ion::math::Dot(ray_end - ray_origin, ray_direction) /
                          ray_direction_length2;
        depths.push_back(new_depth);
      }
    }
  }
  return image::Ldi4f(ldi.GetSize(), sample_counts, std::move(colors),
                      std::move(depths));
}

// Ray-sphere intersection using the ray for the specified |pixel|.
bool IntersectSphere(const Camera& camera, const Point2i& pixel,
                     const SphereViewGroupLoader::Sphere& sphere, float* t_hit,
                     Point3f* hitpoint) {
  Point3f start = camera.RayOrigin(pixel);
  Vector3f direction = camera.RayDirection(pixel);
  if (!geometry::ComputeRaySphereIntersection(sphere.center, sphere.radius,
                                              start, direction, t_hit)) {
    return false;
  }
  *hitpoint = start + direction * (*t_hit);
  return true;
}

// Composites the |rgbd| image onto the given |ldi|, removing occluded
// samples.
//
// The |ldi| and |rgbd| images must have the same depth units.
//
// This will *not* work well if the input Ldi has multiple layers.
Ldi4f CompositeLdiAndRGBD(const Ldi4f& ldi, const Image4f& rgbd) {
  std::vector<int> sample_counts;
  std::vector<Color4f> colors;
  std::vector<float> depths;
  for (int y = 0; y < ldi.GetHeight(); ++y) {
    for (int x = 0; x < ldi.GetWidth(); ++x) {
      CHECK_LE(ldi.GetSampleCount({x, y}), 1)
          << "SphereViewGroupLoader does not (yet) support multi-layer LDIs";
      // Use the samples from the LDI if it has any samples which would
      // occlude the RGBD image.
      //
      // It is sufficient to only check the front-most depth value because
      // LDI samples should be sorted front-to-back.
      bool use_ldi_values = !ldi.GetDepths({x, y}).empty() &&
                            ldi.GetDepths({x, y}).front() < rgbd.At(x, y)[3];
      // Definitely use values from the original LDI if the RGBD image is
      // infinite here.
      use_ldi_values |= !std::isfinite(rgbd.At(x, y)[3]);
      if (use_ldi_values) {
        // Pass through the original sample-counts & colors.
        sample_counts.push_back(ldi.GetSampleCount({x, y}));
        for (int s = 0; s < ldi.GetSampleCount({x, y}); ++s) {
          colors.push_back(ldi.GetColors({x, y})[s]);
          depths.push_back(ldi.GetDepths({x, y})[s]);
        }
      } else {
        // Replace the LDI pixels with the rgbd value.
        Color4f rgbd_value = rgbd.At(x, y);
        Color4f rgba = rgbd_value;
        rgba[3] = 1.0f;
        float depth = rgbd_value[3];

        colors.push_back(rgba);
        depths.push_back(depth);
        sample_counts.push_back(1);
      }
    }
  }
  return image::Ldi4f(ldi.GetSize(), sample_counts, std::move(colors),
                      std::move(depths));
}

}  // namespace

Image4f SphereViewGroupLoader::RenderRGBDForCamera(const Camera& camera) const {
  Image4f image(camera.GetImageSize());
  base::ParallelFor(thread_count_, image.Height(), [&](int y) {
    for (int x = 0; x < image.Width(); ++x) {
      Color4f rgbd(Color3f::Zero(), std::numeric_limits<float>::infinity());
      for (const Sphere& sphere : spheres_) {
        float hit_depth;
        Point3f hit_point;
        if (!IntersectSphere(camera, {x, y}, sphere, &hit_depth, &hit_point)) {
          continue;
        }

        // Take the closest hitpoint.
        if (hit_depth < rgbd[3]) {
          Color3f color = shader_(sphere, hit_point);
          rgbd = Color4f(color, hit_depth);
        }
      }

      image.At(x, y) = rgbd;
    }
  });
  return image;
}

base::Status SphereViewGroupLoader::LoadViewGroup(
    int view_group_index, std::vector<std::shared_ptr<Camera>>* cameras,
    std::vector<Ldi4f>* ldis) const {
  std::vector<std::shared_ptr<base::Camera>> camera_buffer;

  base::Status status =
      original_->LoadViewGroup(view_group_index, &camera_buffer, ldis);

  if (ldis) {
    CHECK_EQ(ldis->size(), camera_buffer.size());
    // Add the spheres here.
    //
    // Note that everything must be converted to RayDepth to ensure we can
    // consistently combine the added spheres to the original LDI.
    for (int view = 0; view < ldis->size(); ++view) {
      Ldi4f ldi = ConvertToRayDepth(ldis->at(view), *camera_buffer[view]);
      Image4f rgbd = RenderRGBDForCamera(RayDepthCamera(camera_buffer[view]));
      ldis->at(view) = CompositeLdiAndRGBD(ldi, rgbd);
    }
  }

  for (auto& camera : camera_buffer) {
    camera = std::make_shared<RayDepthCamera>(camera);
  }

  if (cameras) {
    *cameras = std::move(camera_buffer);
  }

  return status;
}

}  // namespace ingest
}  // namespace seurat
