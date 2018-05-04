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

#include "ion/gfx/shaderinputregistry.h"
#include "ion/gfx/statetable.h"
#include "ion/math/vector.h"

namespace seurat {
namespace viewer {

MonoScene::MonoScene(const ion::math::Range2i& viewport_bounds) {
  ion::gfx::ShaderInputRegistryPtr reg =
      ion::gfx::ShaderInputRegistry::GetGlobalRegistry();
  root_.Reset(new ion::gfx::Node);
  root_->SetLabel("Root");
  root_->SetStateTable(ion::gfx::StateTablePtr(new ion::gfx::StateTable));
  root_->GetStateTable()->SetClearColor(ion::math::Vector4f::Zero());
  root_->GetStateTable()->SetClearDepthValue(1.0f);
  root_->GetStateTable()->SetDepthWriteMask(true);
  root_->GetStateTable()->SetViewport(viewport_bounds);
  root_->AddUniform(reg->Create<ion::gfx::Uniform, int>("uIsRightEye", 0));
}

void MonoScene::AddNode(const ion::gfx::NodePtr& node) {
  root_->AddChild(node);
}

void MonoScene::Clear() { root_->ClearChildren(); }

StereoScene::StereoScene(
    const std::array<ion::math::Range2i, 2>& viewport_bounds) {
  ion::gfx::ShaderInputRegistryPtr reg =
      ion::gfx::ShaderInputRegistry::GetGlobalRegistry();
  root_.Reset(new ion::gfx::Node);
  root_->SetLabel("Root");
  root_->SetStateTable(ion::gfx::StateTablePtr(new ion::gfx::StateTable));
  root_->GetStateTable()->SetClearColor(ion::math::Vector4f::Zero());
  root_->GetStateTable()->SetClearDepthValue(1.0f);
  root_->GetStateTable()->SetDepthWriteMask(true);

  for (int eye = 0; eye < 2; ++eye) {
    eye_roots_[eye].Reset(new ion::gfx::Node);
    eye_roots_[eye]->SetLabel(eye == 0 ? "LeftRoot" : "RightRoot");
    eye_roots_[eye]->SetStateTable(
        ion::gfx::StateTablePtr(new ion::gfx::StateTable));
    eye_roots_[eye]->GetStateTable()->SetViewport(viewport_bounds[eye]);
    eye_roots_[eye]->AddUniform(
        reg->Create<ion::gfx::Uniform, int>("uIsRightEye", eye));
    root_->AddChild(eye_roots_[eye]);
  }
}

void StereoScene::AddNode(const ion::gfx::NodePtr& node) {
  for (int eye = 0; eye < 2; ++eye) {
    eye_roots_[eye]->AddChild(node);
  }
}

void StereoScene::Clear() {
  for (int eye = 0; eye < 2; ++eye) {
    eye_roots_[eye]->ClearChildren();
  }
}

}  // namespace viewer
}  // namespace seurat
