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

#ifndef VR_SEURAT_COMPONENT_SHAPE_COMPONENT_UTILS_H_
#define VR_SEURAT_COMPONENT_SHAPE_COMPONENT_UTILS_H_

#include <memory>
#include <vector>

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "seurat/base/color.h"
#include "seurat/component/component.h"
#include "seurat/component/shape_component.h"

namespace seurat {
namespace component {

// Builds the vertices for a wireframe, axis-aligned box covering the range from
// the |minimal| to the |maximal| point, with given |color|.
std::vector<ShapeComponent::VertexAttributes> MakeWireframeBoxVertices(
    const ion::math::Point3f& minimal, const ion::math::Point3f& maximal,
    const base::Color4f& color);

// Builds a ShapeComponent containing a wirefrime, axis-aligned box covering the
// range from the |minimal| to the |maximal| point, with given |color|.
std::unique_ptr<const ShapeComponent> MakeWireframeBoxShape(
    std::string label, const ion::math::Point3f& minimal,
    const ion::math::Point3f& maximal, const base::Color4f& color);

}  // namespace component
}  // namespace seurat

#endif  // VR_SEURAT_COMPONENT_SHAPE_COMPONENT_UTILS_H_
