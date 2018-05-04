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

#include "seurat/baker/framework/ray_classifier.h"

#include "ion/math/matrix.h"
#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/baker/framework/ray_bundle.h"

namespace seurat {
namespace baker {
namespace {

using base::Color4f;
using geometry::Plane3f;
using ion::math::Matrix4f;
using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Point3f;
using ion::math::Vector2i;
using ion::math::Vector3f;

TEST(RayClassifierTest, ProjectingRayClassifier_CollectSolidSample) {
  ExplicitRayBundle bundle;
  for (int y = 0; y < 3; ++y) {
    for (int x = 0; x < 4; ++x) {
      bundle.AddRay(Point3f(x, y, 0.0f), -Vector3f::AxisZ(),
                    std::array<float, 1>{{10.0f}},
                    std::array<Color4f, 1>{{Color4f()}});
    }
  }

  Frame frame;
  frame.draw_order = 0;
  // A quad at z = -10 which should intersect the ray corresponding to the
  // sample at (2, 1), which is ray #6.
  const Point3f start(1.57f, 0.5f, -10.0f);
  frame.quad = {{start,                              //
                 start + Point3f(1.0f, 0.0f, 0.0f),  //
                 start + Point3f(1.0f, 1.0f, 0.0f),  //
                 start + Point3f(0.0f, 1.0f, 0.0f)}};

  std::vector<Frame> all_frames;
  all_frames.push_back(frame);
  const int kThreadCount = 3;
  const float kSecondaryFrameThreshold = 0.0f;
  ProjectingRayClassifier classifier(
      kThreadCount, ProjectingRayClassifier::RenderingMode::kDrawOrder,
      kSecondaryFrameThreshold);
  classifier.Init(all_frames);

  std::vector<RayClassifier::ClassifiedRays> rays_per_frame =
      classifier.ClassifyRays(bundle);
  EXPECT_EQ(1, rays_per_frame.size());

  const RayClassifier::ClassifiedRays& rays = rays_per_frame.front();
  EXPECT_EQ(1, rays.solid_samples.size());
  EXPECT_EQ(std::make_tuple(6, 0), rays.solid_samples[0]);
  EXPECT_TRUE(rays.freespace_rays.empty());
}

TEST(RayClassifierTest, ProjectingRayClassifier_CollectSampleWithOcclusion) {
  // Construct a view with a single ray from (0, 0, 0) to (0, 0, 10).
  ExplicitRayBundle bundle;
  bundle.AddRay(Point3f(0.0f, 0.0f, 0.0f), Vector3f::AxisZ(),
                std::array<float, 1>{{10.0f}},
                std::array<Color4f, 1>{{Color4f()}});

  Frame primary_frame;
  primary_frame.draw_order = 0;
  // A quad at z = 10 which should intersect the ray.
  const Point3f start(-0.25f, -0.22f, 10.0f);
  primary_frame.quad = {{start,                              //
                         start + Point3f(0.5f, 0.0f, 0.0f),  //
                         start + Point3f(0.5f, 0.5f, 0.0f),  //
                         start + Point3f(0.0f, 0.5f, 0.0f)}};

  const int kThreadCount = 3;
  const float kSecondaryFrameThreshold = 0.5f;
  Frame secondary_frame;
  secondary_frame.draw_order = 1;
  secondary_frame.quad = primary_frame.quad;
  for (auto& pt : secondary_frame.quad) {
    pt[2] += 1.0f;
  }

  Frame occluding_frame;
  occluding_frame.draw_order = 2;
  occluding_frame.quad = primary_frame.quad;
  for (auto& pt : occluding_frame.quad) {
    pt[2] -= 6.0f;
  }

  std::vector<Frame> all_frames;
  all_frames.push_back(primary_frame);
  all_frames.push_back(secondary_frame);
  all_frames.push_back(occluding_frame);

  ProjectingRayClassifier classifier(
      kThreadCount, ProjectingRayClassifier::RenderingMode::kDrawOrder,
      kSecondaryFrameThreshold);
  classifier.Init(all_frames);

  std::vector<RayClassifier::ClassifiedRays> rays_per_frame =
      classifier.ClassifyRays(bundle);
  EXPECT_EQ(3, rays_per_frame.size());

  const RayClassifier::ClassifiedRays& primary_frame_rays = rays_per_frame[0];
  const RayClassifier::ClassifiedRays& secondary_frame_rays = rays_per_frame[1];
  const RayClassifier::ClassifiedRays& occluding_frame_rays = rays_per_frame[2];

  EXPECT_EQ(1, primary_frame_rays.solid_samples.size());
  EXPECT_EQ(0, primary_frame_rays.freespace_rays.size());

  EXPECT_EQ(1, secondary_frame_rays.solid_samples.size());
  EXPECT_EQ(0, secondary_frame_rays.freespace_rays.size());

  EXPECT_EQ(0, occluding_frame_rays.solid_samples.size());
  EXPECT_EQ(1, occluding_frame_rays.freespace_rays.size());
}

TEST(RayClassifierTest, DilatingRayClassifier_CollectSampleWithOcclusion) {
  // Construct a view with a single ray from (-1, -1, 0) to (-1, -1, 10).
  ExplicitRayBundle bundle;
  bundle.AddRay(Point3f(-1.0f, -1.0f, 0.0f), Vector3f::AxisZ(),
                std::array<float, 1>{{10.0f}},
                std::array<Color4f, 1>{{Color4f()}});

  Frame primary_frame;
  primary_frame.draw_order = 0;
  // A quad at z = 10 which should *not* intersect the ray unless dilated.
  //
  // The quad is large enough that quantization to integer-texture-sizes will
  // not be significant.  This only matters for the sake of this test, which
  // uses AreaTextureSizer.
  const Point3f start(-0.25f, -0.22f, 10.0f);
  primary_frame.quad = {{start,                                //
                         start + Point3f(50.0f, 0.0f, 0.0f),   //
                         start + Point3f(50.0f, 50.0f, 0.0f),  //
                         start + Point3f(0.0f, 50.0f, 0.0f)}};

  const int kThreadCount = 3;
  const float kSecondaryFrameThreshold = 0.5f;
  Frame secondary_frame;
  secondary_frame.draw_order = 1;
  secondary_frame.quad = primary_frame.quad;
  for (auto& pt : secondary_frame.quad) {
    pt[2] += 1.0f;
  }

  Frame occluding_frame;
  occluding_frame.draw_order = 2;
  occluding_frame.quad = primary_frame.quad;
  for (auto& pt : occluding_frame.quad) {
    pt[2] -= 6.0f;
  }

  std::vector<Frame> all_frames;
  all_frames.push_back(primary_frame);
  all_frames.push_back(secondary_frame);
  all_frames.push_back(occluding_frame);

  // Construct the DilatingRayClassifier with a filter size of 2.5 and an
  // AreaTextureSizer.
  //
  // This should dilate all of the quads by 2.5 units along the x and y planes,
  // enough to collect the ray in the bundle.
  const float kFilterSize = 2.5f;
  DilatingRayClassifier classifier(
      kFilterSize, std::unique_ptr<TextureSizer>(new AreaTextureSizer),
      std::unique_ptr<RayClassifier>(new ProjectingRayClassifier(
          kThreadCount, ProjectingRayClassifier::RenderingMode::kDrawOrder,
          kSecondaryFrameThreshold)));
  classifier.Init(all_frames);

  std::vector<RayClassifier::ClassifiedRays> rays_per_frame =
      classifier.ClassifyRays(bundle);
  EXPECT_EQ(3, rays_per_frame.size());

  const RayClassifier::ClassifiedRays& primary_frame_rays = rays_per_frame[0];
  const RayClassifier::ClassifiedRays& secondary_frame_rays = rays_per_frame[1];
  const RayClassifier::ClassifiedRays& occluding_frame_rays = rays_per_frame[2];

  EXPECT_EQ(1, primary_frame_rays.solid_samples.size());
  EXPECT_EQ(0, primary_frame_rays.freespace_rays.size());

  EXPECT_EQ(1, secondary_frame_rays.solid_samples.size());
  EXPECT_EQ(0, secondary_frame_rays.freespace_rays.size());

  EXPECT_EQ(0, occluding_frame_rays.solid_samples.size());
  EXPECT_EQ(1, occluding_frame_rays.freespace_rays.size());
}

}  // namespace
}  // namespace baker
}  // namespace seurat
