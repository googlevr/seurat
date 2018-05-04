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

#ifndef VR_SEURAT_ARTIFACT_ARTIFACT_H_
#define VR_SEURAT_ARTIFACT_ARTIFACT_H_

#include <memory>

#include "seurat/component/component.h"
#include "seurat/geometry/mesh.h"
#include "seurat/geometry/quad.h"
#include "seurat/geometry/quad_mesh.h"
#include "seurat/image/image.h"

namespace seurat {
namespace artifact {

// A collection of items which can be exported to disk.
struct Artifact {
  std::shared_ptr<const component::Component> component;
  std::shared_ptr<const image::Image4f> texture;
  std::shared_ptr<const geometry::QuadMesh> quad_mesh;
  std::shared_ptr<const geometry::Mesh> mesh;
};

}  // namespace artifact
}  // namespace seurat

#endif  // VR_SEURAT_ARTIFACT_ARTIFACT_H_
