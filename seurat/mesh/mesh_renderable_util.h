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

#ifndef VR_SEURAT_MESH_MESH_RENDERABLE_UTIL_H_
#define VR_SEURAT_MESH_MESH_RENDERABLE_UTIL_H_

#include <limits>
#include <string>
#include <vector>

#include "ion/base/datacontainer.h"
#include "ion/gfx/bufferobject.h"
#include "ion/gfx/indexbuffer.h"
#include "ion/gfx/shaderinputregistry.h"
#include "ion/gfx/shaderprogram.h"
#include "ion/gfx/shape.h"
#include "ion/gfxutils/shadermanager.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/component/renderable.h"

namespace seurat {
namespace mesh {

// Builds the root node in the scene graph.
ion::gfx::NodePtr BuildRenderNode(
    const component::Renderable::RenderingContext &context,
    const ion::gfx::ShaderInputRegistryPtr &shader_registry,
    const ion::gfx::ShaderProgramPtr &shader_program,
    const ion::gfx::ShapePtr &shape, size_t *u_clip_from_mesh_matrix_index);

// Wraps the given |node| to set up alpha blending.
ion::gfx::NodePtr BuildTransparentPassNode(const ion::gfx::NodePtr &node,
                                           bool is_alpha_premultiplied);

}  // namespace mesh
}  // namespace seurat

#endif  // VR_SEURAT_MESH_MESH_RENDERABLE_UTIL_H_
