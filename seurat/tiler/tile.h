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

#ifndef VR_SEURAT_TILER_TILE_H_
#define VR_SEURAT_TILER_TILE_H_

#include "seurat/geometry/quad.h"

namespace seurat {
namespace tiler {

struct Tile {
  // The Subdivision cell from which this tile was generated.
  //
  // See tiler::Subdivision for details.
  int cell;

  // The geometry of this tile.
  geometry::Quad3f quad;
};

}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_TILE_H_
