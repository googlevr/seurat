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

#ifndef VR_SEURAT_VIEWER_SCENE_H_
#define VR_SEURAT_VIEWER_SCENE_H_

#include <array>

#include "ion/gfx/node.h"
#include "ion/math/range.h"

namespace seurat {
namespace viewer {

// A scene manages the top of the Ion node hierarchy, holding either a single
// root node (MonoScene) or a root and two eye nodes (StereoScene). The scene is
// responsible for setting the viewport and for setting a uniform "uIsRightEye"
// to indicate which eye the sub-tree represents. Mono scenes only have a single
// (left) eye, represented by the root itself.
//
// An example node hierarchy for a MonoScene is shown below. The example nodes
// 'A' and 'B' describe the actual scene content. They would be added by calls
// to AddNode().
//
//      Root
//      /  \
//     A    B
//
// The equivalent hierarchy for a StereoScene has the additional eye nodes 'L'
// and 'R' for the left and right eye, respectively:
//
//      Root
//     /    \
//    L      R
//   / \    / \
//  A   B  A   B

class Scene {
 public:
  virtual ~Scene() = default;
  // Adds a node to the scene.
  virtual void AddNode(const ion::gfx::NodePtr& node) = 0;
  // Returns the root node.
  virtual const ion::gfx::NodePtr& GetRoot() = 0;
  // Removes all previously added nodes.
  virtual void Clear() = 0;
};

class MonoScene : public Scene {
 public:
  explicit MonoScene(const ion::math::Range2i& viewport_bounds);
  void AddNode(const ion::gfx::NodePtr& node) override;
  const ion::gfx::NodePtr& GetRoot() override { return root_; }
  void Clear() override;

 private:
  ion::gfx::NodePtr root_;
};

class StereoScene : public Scene {
 public:
  // Constructs a stereo scene for the given
  // |viewport_bounds|. |viewport_bounds[0]| is the viewport for the left eye
  // and |viewport_bounds[1]| is the viewport for the right eye.
  explicit StereoScene(
      const std::array<ion::math::Range2i, 2>& viewport_bounds);

  void AddNode(const ion::gfx::NodePtr& node) override;
  const ion::gfx::NodePtr& GetRoot() override { return root_; }
  void Clear() override;

 private:
  ion::gfx::NodePtr root_;
  std::array<ion::gfx::NodePtr, 2> eye_roots_;
};

}  // namespace viewer
}  // namespace seurat

#endif  // VR_SEURAT_VIEWER_SCENE_H_
