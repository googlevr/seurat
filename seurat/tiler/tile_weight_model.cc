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

#include "seurat/tiler/tile_weight_model.h"

#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/geometry/bilinear_interpolator.h"
#include "seurat/geometry/fibonacci_sphere.h"
#include "seurat/geometry/plane.h"
#include "seurat/geometry/polygon.h"
#include "seurat/geometry/ray_sphere_intersection.h"
#include "seurat/geometry/triangle.h"
#include "seurat/tiler/tile.h"

namespace seurat {
namespace tiler {

using geometry::BilinearInterpolator;
using geometry::Plane3f;
using geometry::Quad2f;
using geometry::Quad3f;
using geometry::Triangle3f;
using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point3d;
using ion::math::Point3f;
using ion::math::Range2f;
using ion::math::Vector2f;
using ion::math::Vector2i;
using ion::math::Vector3f;

bool ProjectedAreaTileWeightModel::GetTileWeight(
    const Tile& tile, absl::Span<float> weight) const {
  DCHECK_EQ(weight.size(), GetDimension());
  // Project onto the unit-sphere and compute the area of the resulting quad.
  //
  // This assumes that the quad is relatively-small, such that a
  // small-angle approximation holds and the area of the quad is approximately
  // the area of its projection on the unit-sphere.
  Quad3f quad_projected;
  for (int i = 0; i < 4; ++i) {
    quad_projected[i] =
        ion::math::Normalized(tile.quad[i] - Point3f::Zero()) + Point3f::Zero();
  }

  float estimated_area = 0.5f * ion::math::Length(ion::math::Cross(
                                    quad_projected[2] - quad_projected[0],
                                    quad_projected[3] - quad_projected[1]));

  if (!std::isfinite(estimated_area)) {
    return false;
  }
  // Normalize such that 1.0 indicates a full sphere.
  weight[0] = estimated_area / (4.0f * M_PI);
  return true;
}

int CombinedTileWeightModel::GetDimension() const {
  int dim = 0;
  for (const auto& model : models_) {
    dim += model->GetDimension();
  }
  return dim;
}

bool CombinedTileWeightModel::GetTileWeight(const Tile& tile,
                                            absl::Span<float> weight) const {
  DCHECK_EQ(weight.size(), GetDimension());
  int offset = 0;
  for (const auto& model : models_) {
    const int dim = model->GetDimension();
    if (!model->GetTileWeight(tile, weight.subspan(offset, dim))) {
      return false;
    }
    offset += dim;
  }
  return true;
}

std::unique_ptr<TileWeightModel> DirectionalOverdrawTileWeightModel::Build(
    int samples, float field_of_view, float headbox_radius) {
  std::vector<Vector3f> directions;
  directions.reserve(samples);
  for (int s = 0; s < samples; ++s) {
    directions.push_back(
        Vector3f(geometry::GenerateFibonacciSpherePoint(samples, 0.0, s) -
                 Point3d::Zero()));
  }
  return std::unique_ptr<TileWeightModel>(
      new DirectionalOverdrawTileWeightModel(std::move(directions),
                                             field_of_view, headbox_radius));
}

namespace {

// Determines how much to subdivide the given |quad| in each direction to reach
// the specified number of subdivisions per revolution.
Vector2i GetSubdivisionFactor(const Quad3f& quad,
                              int subdivisions_per_revolution) {
  // Compute the angle formed by vertices of the quad & the origin.
  float angle_x_radians =
      std::max(std::acos(ion::math::Dot(
                   ion::math::Normalized(quad[0] - Point3f::Zero()),
                   ion::math::Normalized(quad[1] - Point3f::Zero()))),
               std::acos(ion::math::Dot(
                   ion::math::Normalized(quad[2] - Point3f::Zero()),
                   ion::math::Normalized(quad[3] - Point3f::Zero()))));
  float angle_y_radians =
      std::max(std::acos(ion::math::Dot(
                   ion::math::Normalized(quad[0] - Point3f::Zero()),
                   ion::math::Normalized(quad[3] - Point3f::Zero()))),
               std::acos(ion::math::Dot(
                   ion::math::Normalized(quad[1] - Point3f::Zero()),
                   ion::math::Normalized(quad[2] - Point3f::Zero()))));
  Vector2i subdivision_resolution(Vector2f(angle_x_radians, angle_y_radians) *
                                  subdivisions_per_revolution /
                                  (2.0f * static_cast<float>(M_PI)));
  subdivision_resolution[0] = std::max(1, subdivision_resolution[0]);
  subdivision_resolution[1] = std::max(1, subdivision_resolution[1]);
  return subdivision_resolution;
}

// Returns the closest point within an origin-centered sphere of the specified
// radius to the ray.
Point3f ClosestPointInSphere(const Point3f& ray_start,
                             const Vector3f& ray_direction_normalized,
                             float radius) {
  float t_hit;
  float squared_radius = radius * radius;
  if (ion::math::LengthSquared(ray_start - Point3f::Zero()) < squared_radius) {
    // If the ray starts inside the sphere, return the starting point.
    return ray_start;
  }
  if (geometry::ComputeRaySphereIntersection(Point3f::Zero(), radius, ray_start,
                                             ray_direction_normalized,
                                             &t_hit)) {
    return ray_start + ray_direction_normalized * t_hit;
  } else {
    // If the ray misses the sphere, find the closest point on the sphere.

    // Project the vector from ray_start->origin onto ray_direction.
    Vector3f projection =
        ion::math::Dot(Point3f::Zero() - ray_start, ray_direction_normalized) *
        ray_direction_normalized;

    // Add this to the ray_start to get the point along the ray which is
    // closest to the origin.
    Point3f closest_point_to_origin = projection + ray_start;

    // Scale this point to get the closest point on the sphere.
    return ion::math::Normalized(closest_point_to_origin - Point3f::Zero()) *
               radius +
           Point3f::Zero();
  }
}

}  // namespace

bool DirectionalOverdrawTileWeightModel::GetTileWeight(
    const Tile& tile, absl::Span<float> weight) const {
  // The tile weight is computed as follows:
  //
  //  * Subdivide the tile's quad into tiny pieces, according to a target
  //    subdivisions-per-revolution.  This improves accuracy of the weight,
  //    which is essentially sampling an integral over the quad.
  //  * For each patch, find the eye point (within the viewing volume), which
  //    has the "worst case" view of the quad:  the most head-on point-of-view.
  //  * Compute the projected-area of that patch.
  //  * If the patch is seen by a cone-shaped frustum facing a particular
  //     direction, then use this projected-area as the weight (normalized
  //     according the field-of-view of the cone).

  DCHECK_EQ(weight.size(), GetDimension());

  // This computes the surface area of the spherical cap of the circular field
  // of view used to model directional overdraw.
  //
  // Based on Wikipedia/Steradian.
  const float steradians = 2.0f * M_PI * (1.0f - cos_half_field_of_view_);

  std::fill(weight.begin(), weight.end(), 0.0f);

  // Subdivide into bilinearly-interpolated patches for integration.
  const int kSubdivisionsPerRevolution = 100;

  Vector2i subdivision_factor =
      GetSubdivisionFactor(tile.quad, kSubdivisionsPerRevolution);
  Point2f subdivision_scale =
      Point2f(1.0f / subdivision_factor[0], 1.0f / subdivision_factor[1]);

  const Vector3f normal = geometry::NormalFromTriangle(
      Triangle3f{{tile.quad[0], tile.quad[1], tile.quad[2]}});

  BilinearInterpolator<Point3f> tile_interp(tile.quad);
  for (int x = 0; x < subdivision_factor[0]; ++x) {
    for (int y = 0; y < subdivision_factor[1]; ++y) {
      const Quad3f patch = {
          {tile_interp.At(Point2f(x, y) * subdivision_scale),
           tile_interp.At(Point2f(x + 1, y) * subdivision_scale),
           tile_interp.At(Point2f(x + 1, y + 1) * subdivision_scale),
           tile_interp.At(Point2f(x, y + 1) * subdivision_scale)}};

      Point3f average =
          tile_interp.At(Point2f::Zero() + (Vector2f(x + 0.5f, y + 0.5f) /
                                            Vector2f(subdivision_factor)));
      Vector3f average_dir = ion::math::Normalized(average - Point3f::Zero());

      // Select the worst-case eye position (based on headbox size) which
      // maximizes projected area.
      const Point3f eye =
          ClosestPointInSphere(average, normal, headbox_radius_);
      std::array<Point3f, 4> quad_projected;
      for (int i = 0; i < 4; ++i) {
        quad_projected[i] =
            ion::math::Normalized(patch[i] - eye) + Point3f::Zero();
      }
      // Area estimated with the small-angle approximation, which should hold
      // since this is subdivided specifically to encourage small angles here.
      float estimated_area = 0.5f * ion::math::Length(ion::math::Cross(
                                        quad_projected[2] - quad_projected[0],
                                        quad_projected[3] - quad_projected[1]));
      float estimated_area_normalized = estimated_area / steradians;

      for (int s = 0; s < directions_.size(); ++s) {
        // Only add this weight if the center of the patch is within the cone
        // centered around directions_[s] with the specified field-of-view.
        if (ion::math::Dot(average_dir, directions_[s]) >
            cos_half_field_of_view_) {
          weight[s] += estimated_area_normalized;
        }
      }
    }
  }

  return true;
}

}  // namespace tiler
}  // namespace seurat
