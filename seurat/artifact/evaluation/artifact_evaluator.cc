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

#include "seurat/artifact/evaluation/artifact_evaluator.h"

#include <mutex>  // NOLINT(build/c++11)
#include <string>
#include <vector>

#include "ion/math/vector.h"
#include "absl/strings/str_cat.h"
#include "seurat/base/color.h"
#include "seurat/base/parallel.h"
#include "seurat/base/reporting.h"
#include "seurat/geometry/mesh.h"
#include "seurat/geometry/raytracer.h"
#include "seurat/image/image.h"
#include "seurat/image/ldi.h"
#include "seurat/image/ldi_util.h"
#include "seurat/ingest/view_group_loader_util.h"

namespace seurat {
namespace artifact {

using base::Color3f;
using base::Color4f;
using geometry::Mesh;
using geometry::Raytracer;
using image::Image3f;
using image::Image4f;
using image::Ldi4f;
using ion::math::Point3f;
using ion::math::Vector2i;
using ion::math::Vector3f;

base::Status EvaluationExporter::Process(Artifact* artifact) const {
  // Print the results of the evaluation in the following format:
  // [Evaluator Name]
  //   [Result Name]:  [Result Value]
  //   ...
  // ...

  const char* const kIndent = " ";
  std::string results_string;
  for (const auto& evaluator : evaluators_) {
    results_string += evaluator->Name();
    results_string += "\n";
    std::vector<ArtifactEvaluator::Result> results =
        evaluator->Evaluate(*artifact);
    for (const auto& result : results) {
      results_string += kIndent;
      results_string += result.name;
      results_string += ": ";
      results_string += result.value;
      results_string += "\n";
    }
  }
  std::string full_filename = absl::StrCat(basename_, ".eval");
  return filesystem_->SetContents(full_filename, results_string);
}

std::vector<ArtifactEvaluator::Result> OverdrawEvaluator::Evaluate(
    const Artifact& artifact) const {
  const char* const kMinOverdraw = "single view min overdraw";
  const char* const kMaxOverdraw = "single view max overdraw";
  const char* const kAvgOverdraw = "all view avg overdraw";

  std::vector<Result> results;

  std::vector<CostEstimator::SceneGeometry> geometry;
  if (!artifact.mesh) {
    return results;
  }
  const Mesh& mesh = *artifact.mesh;
  geometry.push_back({mesh, 1.0f, false /* opaque=false */});

  CostEstimator cost_estimator = CostEstimator::Build(geometry);

  std::vector<std::shared_ptr<base::Camera>> all_cameras;
  for (int i = 0; i < view_loader_->GetNumViewGroups(); ++i) {
    std::vector<std::shared_ptr<base::Camera>> cameras;
    base::Status result = view_loader_->LoadViewGroup(i, &cameras, nullptr);
    if (!result.ok()) {
      base::SeuratError(result.error_message());
      continue;
    }
    all_cameras.insert(all_cameras.end(), cameras.begin(), cameras.end());
  }

  std::vector<float> overdraw_factor_per_camera(all_cameras.size());

  base::ParallelFor(thread_count_, all_cameras.size(), [&](int camera_index) {
    const std::shared_ptr<base::Camera>& camera = all_cameras[camera_index];
    Vector2i size = camera->GetImageSize();
    float total_cost = 0.0f;
    for (int y = 0; y < size[1]; ++y) {
      for (int x = 0; x < size[0]; ++x) {
        auto origin = camera->RayOrigin({x, y});
        auto direction = camera->RayDirection({x, y});
        total_cost += cost_estimator.EstimateRayCost(origin, direction);
      }
    }
    float avg_cost = total_cost / (size[0] * size[1]);
    overdraw_factor_per_camera[camera_index] = avg_cost;
  });

  float min_overdraw = *std::min_element(overdraw_factor_per_camera.begin(),
                                         overdraw_factor_per_camera.end());
  float max_overdraw = *std::max_element(overdraw_factor_per_camera.begin(),
                                         overdraw_factor_per_camera.end());
  float avg_overdraw = std::accumulate(overdraw_factor_per_camera.begin(),
                                       overdraw_factor_per_camera.end(), 0.0f) /
                       static_cast<float>(overdraw_factor_per_camera.size());

  results.push_back({kMinOverdraw, std::to_string(min_overdraw)});
  results.push_back({kMaxOverdraw, std::to_string(max_overdraw)});
  results.push_back({kAvgOverdraw, std::to_string(avg_overdraw)});

  return results;
}

std::vector<ArtifactEvaluator::Result> RenderEvaluator::Evaluate(
    const Artifact& artifact) const {
  const char* const kMinMse = "single view min MSE";
  const char* const kMaxMse = "single view max MSE";
  const char* const kGlobalMse = "all view MSE";
  const char* const kGlobalPsnr = "all view PSNR";
  const char* const kNonOpaquePixels = "non-opaque pixels (\"cracks\")";
  const char* const kTextureSize = "texture size (pixels)";
  const char* const kTextureOccupancy =
      "texture occupancy (fraction of texels >0.5 alpha)";
  const char* const kMeshSize = "mesh size (triangles)";

  std::vector<Result> results;

  if (!artifact.mesh || !artifact.texture) {
    return results;
  }
  const Mesh& mesh = *artifact.mesh;
  const Image4f& texture = *artifact.texture;

  const int total_pixels = texture.Width() * texture.Height();
  const int total_triangles = mesh.GetTriangleCount();
  int opaque_pixels = 0;
  for (const auto& rgba : texture) {
    if (rgba[3] > 0.5) opaque_pixels++;
  }
  double texture_occupancy = static_cast<double>(opaque_pixels) / total_pixels;
  results.push_back({kTextureSize, std::to_string(total_pixels)});
  results.push_back({kTextureOccupancy, std::to_string(texture_occupancy)});
  results.push_back({kMeshSize, std::to_string(total_triangles)});

  render_sim_->Build(mesh, &texture);

  std::vector<float> mse_per_view;
  std::vector<int> non_opaque_pixels_per_view;

  // Scratch space for each rendered image.
  Image4f simulated_render;

  base::Status status = ingest::ForEachViewGroupPrefetching(
      *view_loader_, [&](std::vector<std::shared_ptr<base::Camera>> cameras,
                         std::vector<Ldi4f> ldis) {
        const int num_views = ldis.size();
        CHECK_EQ(cameras.size(), num_views);

        for (int view = 0; view < num_views; ++view) {
          const std::shared_ptr<base::Camera>& camera = cameras[view];
          Vector2i size = camera->GetImageSize();
          if (size[0] * size[1] == 0) {
            // Skip this view.
            continue;
          }

          render_sim_->Render(*camera, &simulated_render);
          CHECK_EQ(size, simulated_render.GetSize());

          // Compute the flattened ldi to compare with.
          const Image3f expected_image = image::FlattenLdi(ldis[view]);
          CHECK_EQ(size, expected_image.GetSize());

          // Sum of squared differences.
          float ssd = 0.0f;
          int non_opaque_pixels = 0;
          for (int y = 0; y < size[1]; ++y) {
            for (int x = 0; x < size[0]; ++x) {
              const float kMinOpaqueAlpha = 1.0f - (0.5f / 255.0f);
              if (simulated_render.At(x, y)[3] < kMinOpaqueAlpha) {
                non_opaque_pixels++;
                continue;
              }
              Color3f simulated_color(simulated_render.At(x, y)[0],
                                      simulated_render.At(x, y)[1],
                                      simulated_render.At(x, y)[2]);
              Color3f expected_color = expected_image.At(x, y);
              Color3f diff = simulated_color - expected_color;
              float pixel_ssd = 0.0f;
              for (int c = 0; c < 3; ++c) {
                pixel_ssd += diff[c] * diff[c];
              }
              ssd += pixel_ssd / 3.0f;
            }
          }
          int opaque_pixels = size[0] * size[1] - non_opaque_pixels;
          // Verify no divide by zero.
          if (opaque_pixels > 0) {
            float mse = ssd / static_cast<float>(opaque_pixels);
            CHECK(std::isfinite(mse));
            mse_per_view.push_back(mse);
          }
          non_opaque_pixels_per_view.push_back(non_opaque_pixels);
        }
        return base::OkStatus();
      });
  if (!status.ok()) {
    base::SeuratError(status.error_message());
  }

  if (!mse_per_view.empty()) {
    float min_mse = *std::min_element(mse_per_view.begin(), mse_per_view.end());
    float max_mse = *std::max_element(mse_per_view.begin(), mse_per_view.end());
    float global_mse =
        std::accumulate(mse_per_view.begin(), mse_per_view.end(), 0.0) /
        static_cast<double>(mse_per_view.size());
    float global_psnr =
        10.0f * std::log10(255.0f * 255.0f / (global_mse * 255.0f));

    results.push_back({kMinMse, std::to_string(min_mse)});
    results.push_back({kMaxMse, std::to_string(max_mse)});
    results.push_back({kGlobalMse, std::to_string(global_mse)});
    results.push_back({kGlobalPsnr, std::to_string(global_psnr)});
  }
  if (!non_opaque_pixels_per_view.empty()) {
    int total_non_opaque = std::accumulate(non_opaque_pixels_per_view.begin(),
                                           non_opaque_pixels_per_view.end(), 0);
    results.push_back({kNonOpaquePixels, std::to_string(total_non_opaque)});
  }
  return results;
}

std::vector<ArtifactEvaluator::Result> GeometryDistortionEvaluator::Evaluate(
    const Artifact& artifact) const {
  const char* const kMSEWarp = "mean squared sample distortion";
  const char* const kDiscardedPixels = "percent pixels discarded";

  // The scale factor applied to distortion.
  //
  // This is somewhat arbitrary, but is required to ensure the final numbers are
  // large enough to be easily interpreted/compared in the human-readable
  // output.
  const float kScale = 2.0 * M_PI * 10000.0;

  std::vector<Result> results;
  if (!artifact.mesh) {
    return results;
  }
  const Mesh& mesh = *artifact.mesh;

  std::unique_ptr<Raytracer> raytracer = Raytracer::Build(mesh);

  // The total screen-space warping of all pixels.
  //
  // This is measured as the sum of squares of the great-circle-distance between
  // original samples & warped samples, relative to a circle centered on the
  // original camera which observed the sample.
  //
  // Samples are "warped" by projection onto the mesh in the direction of the
  // origin.
  //
  // Values are scaled up (somewhat arbitrarily) by kScale.
  float total_squared_pixel_warp = 0.0f;
  // The number of pixels which could not be projected onto the mesh (i.e. the
  // mesh has holes!).
  int total_discarded = 0;
  // The number of pixels.
  int total_pixels = 0;
  // Mutex for total_pixel_warp, total_pixel_count, and total_pixels_discarded.
  std::mutex total_mutex;

  for (int i = 0; i < view_loader_->GetNumViewGroups(); ++i) {
    std::vector<Ldi4f> ldis;
    std::vector<std::shared_ptr<base::Camera>> cameras;
    base::Status result = view_loader_->LoadViewGroup(i, &cameras, &ldis);
    if (!result.ok()) {
      base::SeuratError(result.error_message());
      continue;
    }

    const int num_views = ldis.size();
    CHECK_EQ(cameras.size(), num_views);

    for (int view = 0; view < num_views; ++view) {
      const std::shared_ptr<base::Camera>& camera = cameras[view];
      Vector2i size = camera->GetImageSize();
      if (size[0] * size[1] == 0) {
        // Skip this view.
        continue;
      }

      const Point3f eye_point =
          ion::math::ProjectPoint(camera->GetWorldFromEye(), Point3f::Zero());

      // Loop over all samples of the LDI, parallelizing over y.
      base::BalancedParallelFor(thread_count_, size[1], [&](int y) {
        std::vector<Raytracer::Intersection> intersections;
        float squared_pixel_warp_scanline = 0.0f;
        int num_discarded_scanline = 0;
        for (int x = 0; x < size[0]; ++x) {
          for (float depth : ldis[view].GetDepths({x, y})) {
            // For this sample (given by {x, y, depth}), compute distortion by
            //  1. Projecting the point onto the mesh (in the direction of the
            //     origin)
            //  2. Computing great-circle-distance between the original point &
            //     projected point, relative to 'eye_point'.

            // depth in the ldi is not necessarily in world-space.
            Point3f end = camera->RayEnd({x, y}, depth);

            // Trace a ray from the origin->end.
            //
            // The resulting intersection which minimizes abs(t_hit - 1.0) is
            // closest to 'end'.
            raytracer->FindAllIntersections(
                Point3f::Zero(), end - Point3f::Zero(), &intersections);

            if (intersections.empty()) {
              num_discarded_scanline++;
              continue;
            }

            // Find the intersection which was closest to the original end
            // point.
            int best_intersection_index = 0;
            for (int intersection_index = 1;
                 intersection_index < intersections.size();
                 ++intersection_index) {
              if (std::fabs(intersections[intersection_index].t_hit - 1.0f) <
                  std::fabs(intersections[best_intersection_index].t_hit -
                            1.0f)) {
                best_intersection_index = intersection_index;
              }
            }

            Point3f warped_end =
                Point3f::Zero() +
                (end - Point3f::Zero()) *
                    intersections[best_intersection_index].t_hit;

            // Distance between 'end' and 'warped_end' with respect to
            // eye-point.
            float great_circle_distance = std::acos(
                ion::math::Dot(ion::math::Normalized(warped_end - eye_point),
                               ion::math::Normalized(end - eye_point)));
            if (!std::isfinite(great_circle_distance)) {
              num_discarded_scanline++;
              continue;
            }

            // Scale up to get values which are easier to interpret in the final
            // human-readable report.
            great_circle_distance *= kScale;

            squared_pixel_warp_scanline +=
                great_circle_distance * great_circle_distance;
          }
        }
        std::lock_guard<std::mutex> lock(total_mutex);
        total_squared_pixel_warp += squared_pixel_warp_scanline;
        total_discarded += num_discarded_scanline;
        total_pixels += size[0];
      });
    }
  }

  float mse_warp = total_squared_pixel_warp / (total_pixels - total_discarded);
  results.push_back(Result{kMSEWarp, std::to_string(mse_warp)});

  float percent_discarded = 100.0f * total_discarded / total_pixels;
  results.push_back(
      Result{kDiscardedPixels, std::to_string(percent_discarded)});

  return results;
}

}  // namespace artifact
}  // namespace seurat
