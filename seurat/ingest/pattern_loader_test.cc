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

#include "seurat/ingest/pattern_loader.h"

#include "ion/math/vector.h"
#include "gtest/gtest.h"

namespace seurat {
namespace ingest {

namespace {

constexpr int kImageSize = 20;
constexpr int kFeatureSize = 4;
constexpr float kAngle = 0.0f;
constexpr float kDepth = 0.2f;

std::unique_ptr<ViewGroupLoader> BuildPatternLoader(
    const std::string& pattern) {
  PatternLoader::Parameters parameters{pattern, kImageSize, kFeatureSize,
                                       kAngle, kDepth};
  return std::unique_ptr<ViewGroupLoader>(new PatternLoader(parameters));
}

TEST(PatternLoader, GetNumViewGroups) {
  const std::array<std::string, 2> patterns = {{"checkerboard", "stripes"}};
  for (const std::string& pattern : patterns) {
    auto view_group_loader = BuildPatternLoader(pattern);
    EXPECT_EQ(1, view_group_loader->GetNumViewGroups());
  }
}

TEST(PatternLoader, LoadViewGroup) {
  const std::array<std::string, 2> patterns = {{"checkerboard", "stripes"}};
  for (const std::string& pattern : patterns) {
    auto view_group_loader = BuildPatternLoader(pattern);
    std::vector<std::shared_ptr<base::Camera>> cameras;
    std::vector<image::Ldi4f> ldis;
    for (int i = 0; i < view_group_loader->GetNumViewGroups(); ++i) {
      EXPECT_TRUE(view_group_loader->LoadViewGroup(i, &cameras, &ldis).ok());
    }
    EXPECT_EQ(1, cameras.size());
    EXPECT_EQ(1, ldis.size());
    EXPECT_EQ(kImageSize, ldis[0].GetWidth());
    EXPECT_EQ(kImageSize, ldis[0].GetHeight());
    EXPECT_EQ(kImageSize * kImageSize, ldis[0].GetSampleCount());

    auto depths = ldis[0].GetDepths();
    for (auto depth : depths) {
      EXPECT_NEAR(kDepth, depth, 1.0e-5f);
    }
  }
}

}  // namespace

}  // namespace ingest
}  // namespace seurat
