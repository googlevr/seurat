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

#include "seurat/viewer/scene.h"

#include <array>

#include "ion/gfx/node.h"
#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "gtest/gtest.h"

namespace seurat {
namespace viewer {
namespace {

using ion::math::Point2i;
using ion::math::Range2i;

TEST(MonoScene, Construction) {
  const Range2i kViewportBounds(Point2i(0, 0), Point2i(200, 100));
  MonoScene scene(kViewportBounds);
  EXPECT_EQ(std::string("Root"), scene.GetRoot()->GetLabel());
  // Mono scene root node is expected to have no children.
  EXPECT_EQ(0, scene.GetRoot()->GetChildren().size());
}

TEST(MonoScene, AddNodeAndClear) {
  const Range2i kViewportBounds(Point2i(0, 0), Point2i(200, 100));
  MonoScene scene(kViewportBounds);
  ion::gfx::NodePtr node_a(new ion::gfx::Node);
  node_a->SetLabel("A");
  scene.AddNode(node_a);
  EXPECT_EQ(1, scene.GetRoot()->GetChildren().size());
  EXPECT_EQ(std::string("A"), scene.GetRoot()->GetChildren()[0]->GetLabel());
  scene.Clear();
  EXPECT_EQ(0, scene.GetRoot()->GetChildren().size());
}

TEST(StereoScene, Construction) {
  const std::array<Range2i, 2> kViewportBounds{
      {Range2i(Point2i(0, 0), Point2i(100, 100)),
       Range2i(Point2i(100, 0), Point2i(200, 100))}};
  StereoScene scene(kViewportBounds);
  EXPECT_EQ(std::string("Root"), scene.GetRoot()->GetLabel());
  // Stereoscene root node is expected to have two children.
  EXPECT_EQ(2, scene.GetRoot()->GetChildren().size());
}

TEST(StereoScene, AddNodeAndClear) {
  const std::array<Range2i, 2> kViewportBounds{
      {Range2i(Point2i(0, 0), Point2i(100, 100)),
       Range2i(Point2i(100, 0), Point2i(200, 100))}};
  StereoScene scene(kViewportBounds);
  ion::gfx::NodePtr node_a(new ion::gfx::Node);
  node_a->SetLabel("A");
  scene.AddNode(node_a);
  EXPECT_EQ(2, scene.GetRoot()->GetChildren().size());
  ion::gfx::NodePtr left = scene.GetRoot()->GetChildren()[0];
  ion::gfx::NodePtr right = scene.GetRoot()->GetChildren()[1];
  EXPECT_EQ(1, left->GetChildren().size());
  EXPECT_EQ(1, right->GetChildren().size());
  EXPECT_EQ(std::string("A"), left->GetChildren()[0]->GetLabel());
  EXPECT_EQ(std::string("A"), right->GetChildren()[0]->GetLabel());
  scene.Clear();
  EXPECT_EQ(0, left->GetChildren().size());
  EXPECT_EQ(0, right->GetChildren().size());
  // Check that eye roots are not deleted.
  EXPECT_EQ(2, scene.GetRoot()->GetChildren().size());
}

}  // namespace
}  // namespace viewer
}  // namespace seurat
