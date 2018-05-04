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

#include "seurat/mesh/mesh_component.h"

#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/component/component_tests_utils.h"

using ion::math::Point3f;

namespace seurat {
namespace mesh {
namespace {

TEST(MeshComponentTest, TestIo) {
  ion::gfx::ImagePtr texture_a =
      base::CreateImage(ion::gfx::Image::kRgba8ui, {16, 16});

  ion::gfx::ImagePtr texture_b =
      base::CreateImage(ion::gfx::Image::kRgba8ui, {8, 8});

  const auto mesh_a = MeshComponent::Create(
      "mesh_a", texture_a,
      {
          {Point3f(0.0f, 0.0f, 0.0f), Point3f(0.0f, 0.0f, 1.0f)},
          {Point3f(1.0f, 0.0f, 0.0f), Point3f(1.0f, 0.0f, 1.0f)},
          {Point3f(1.0f, 1.0f, 0.0f), Point3f(1.0f, 1.0f, 1.0f)},
      },
      {0, 1, 2});
  const auto mesh_b = MeshComponent::Create(
      "mesh_b", texture_b,
      {
          {Point3f(0.0f, 0.0f, 2.0f), Point3f(0.0f, 0.0f, 1.0f)},
          {Point3f(1.0f, 0.0f, 2.0f), Point3f(1.0f, 0.0f, 1.0f)},
          {Point3f(1.0f, 1.0f, 2.0f), Point3f(1.0f, 1.0f, 1.0f)},
      },
      {0, 1, 2});

  component::ComponentCommonTests::TestIo(*mesh_a, *mesh_b);
}

}  // namespace
}  // namespace mesh
}  // namespace seurat
