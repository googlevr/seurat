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

#include "seurat/baker/framework/silhouette.h"

#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace seurat {
namespace baker {
namespace {

using ion::math::Point2f;

TEST(ImplicitSilhouetteTest, TestEmpty) {
  ImplicitSilhouette silhouette({}, {});
  // This could probably go either way, but consistency is good.
  EXPECT_FALSE(silhouette.IsSolidAtPoint(Point2f::Zero()));
}

TEST(ImplicitSilhouetteTest, TestSingleSolidSample) {
  ImplicitSilhouette silhouette({Point2f::Zero()}, {});
  EXPECT_TRUE(silhouette.IsSolidAtPoint(Point2f(1.0f, 1.0f)));
}

TEST(ImplicitSilhouetteTest, TestSingleFreespaceSample) {
  ImplicitSilhouette silhouette({}, {Point2f::Zero()});
  EXPECT_FALSE(silhouette.IsSolidAtPoint(Point2f(1.0f, 1.0f)));
}

TEST(ImplicitSilhouetteTest, TestMultipleSamples) {
  Point2f solid_sample(1.0f, 1.0f);
  Point2f freespace_sample(2.0f, 2.0f);
  ImplicitSilhouette silhouette({solid_sample}, {freespace_sample});

  EXPECT_TRUE(silhouette.IsSolidAtPoint(Point2f(0.1f, 0.1f)));
  EXPECT_TRUE(silhouette.IsSolidAtPoint(Point2f(1.0f, 1.0f)));
  EXPECT_TRUE(silhouette.IsSolidAtPoint(Point2f(1.4f, 1.4f)));
  EXPECT_FALSE(silhouette.IsSolidAtPoint(Point2f(1.6f, 1.6f)));
  EXPECT_FALSE(silhouette.IsSolidAtPoint(Point2f(2.0f, 2.0f)));
  EXPECT_FALSE(silhouette.IsSolidAtPoint(Point2f(2.1f, 2.1f)));
}

TEST(SimpleSilhouetteBufferTest, TestMultipleSamples) {
  Point2f solid_sample(1.0f, 1.0f);
  Point2f freespace_sample(2.0f, 2.0f);

  SimpleSilhouetteBuffer buffer;
  buffer.AddSolidSample(solid_sample);
  buffer.AddFreespaceSample(freespace_sample);

  std::unique_ptr<ImplicitSilhouette> silhouette = buffer.Resolve();

  EXPECT_TRUE(silhouette->IsSolidAtPoint(Point2f(0.1f, 0.1f)));
  EXPECT_TRUE(silhouette->IsSolidAtPoint(Point2f(1.0f, 1.0f)));
  EXPECT_TRUE(silhouette->IsSolidAtPoint(Point2f(1.4f, 1.4f)));
  EXPECT_FALSE(silhouette->IsSolidAtPoint(Point2f(1.6f, 1.6f)));
  EXPECT_FALSE(silhouette->IsSolidAtPoint(Point2f(2.0f, 2.0f)));
  EXPECT_FALSE(silhouette->IsSolidAtPoint(Point2f(2.1f, 2.1f)));
}

TEST(CompactSilhouetteBufferTest, TestMultipleSamples) {
  const float kEpsilon = 1e-6f;

  CompactSilhouetteBuffer buffer({10, 10});

  // These values should not change as a result of quantization to the 10x10
  // bins over [0, 1]^2.
  Point2f solid_sample(0.05f, 0.05f);
  Point2f freespace_sample(0.95f, 0.95f);

  buffer.AddSolidSample(solid_sample);
  buffer.AddFreespaceSample(freespace_sample);

  std::unique_ptr<ImplicitSilhouette> silhouette = buffer.Resolve();

  EXPECT_TRUE(silhouette->IsSolidAtPoint(Point2f(0.0f, 0.0f)));
  EXPECT_TRUE(silhouette->IsSolidAtPoint(Point2f(0.1f, 0.1f)));
  EXPECT_TRUE(
      silhouette->IsSolidAtPoint(Point2f(0.5f - kEpsilon, 0.5f - kEpsilon)));
  EXPECT_FALSE(silhouette->IsSolidAtPoint(Point2f(1.0f, 1.0f)));
  EXPECT_FALSE(silhouette->IsSolidAtPoint(Point2f(0.8f, 0.8f)));
  EXPECT_FALSE(
      silhouette->IsSolidAtPoint(Point2f(0.5f + kEpsilon, 0.5f + kEpsilon)));
}

}  // namespace
}  // namespace baker
}  // namespace seurat
