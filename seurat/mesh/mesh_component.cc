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

#include <array>

#include "ion/math/matrix.h"
#include "ion/math/matrixutils.h"
#include "ion/math/transformutils.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/mesh/mesh_renderable.h"

namespace seurat {
namespace mesh {

using base::StructureSink;
using base::StructureSource;
using component::Component;
using component::Renderable;
using ion::math::Matrix3f;
using ion::math::Point2f;
using ion::math::Point3f;
using ion::math::Range3f;
using ion::math::Vector2f;
using ion::math::Vector3f;

namespace {

Range3f ComputeBoundingBox(
    const std::vector<MeshComponent::VertexAttributes>& verts) {
  Range3f aabb;
  for (const auto& vert : verts) {
    aabb.ExtendByPoint(vert.position);
  }
  return aabb;
}

}  // namespace

// static
std::unique_ptr<MeshComponent> MeshComponent::Create(
    std::string label, ion::gfx::ImagePtr texture_atlas,
    std::vector<VertexAttributes> attribute_buf,
    std::vector<uint32> index_buf) {
  const auto aabb = ComputeBoundingBox(attribute_buf);
  return std::unique_ptr<MeshComponent>(
      new MeshComponent(std::move(label), texture_atlas,
                        std::move(attribute_buf), std::move(index_buf), aabb));
}

// static
std::unique_ptr<const Component> MeshComponent::Create(
    std::string label, StructureSource* source) {
  const int8 has_texture = source->ReadPod<int8>();
  ion::gfx::ImagePtr texture(nullptr);
  if (has_texture) {
    texture = source->ReadImage();
  }

  std::vector<VertexAttributes> attribute_buf(source->ReadPod<int64>());
  for (VertexAttributes& attr : attribute_buf) {
    attr.position = source->ReadPoint<Point3f>();
    attr.tex_coord = source->ReadPoint<Point3f>();
  }

  std::vector<uint32> index_buf(source->ReadPod<int64>());
  source->ReadPodArray(index_buf.data(), index_buf.size());

  const auto aabb_min = source->ReadPoint<Point3f>();
  const auto aabb_max = source->ReadPoint<Point3f>();
  const Range3f aabb(aabb_min, aabb_max);

  return std::unique_ptr<MeshComponent>(
      new MeshComponent(std::move(label), texture, std::move(attribute_buf),
                        std::move(index_buf), aabb));
}

Renderable* MeshComponent::GetRenderable() const { return renderable_.get(); }

bool MeshComponent::operator==(const Component& other) const {
  if (!IsComponentEqualTo(other)) {
    return false;
  }
  const MeshComponent* other_mesh = static_cast<const MeshComponent*>(&other);
  if (other_mesh->GetIndexBuffer() != GetIndexBuffer()) {
    return false;
  }
  if (other_mesh->GetAttributeBuffer() != GetAttributeBuffer()) {
    return false;
  }
  if (other_mesh->GetBoundingBox() != GetBoundingBox()) {
    return false;
  }
  return base::CompareImagesEqual(GetTextureAtlas(),
                                  other_mesh->GetTextureAtlas());
}

MeshComponent::MeshComponent(std::string label,
                             ion::gfx::ImagePtr texture_atlas,
                             std::vector<VertexAttributes> attribute_buf,
                             std::vector<uint32> index_buf,
                             const ion::math::Range3f& aabb)
    : Component(std::move(label)),
      texture_atlas_(texture_atlas),
      attribute_buf_(std::move(attribute_buf)),
      index_buf_(std::move(index_buf)),
      aabb_(aabb),
      renderable_(new MeshRenderable(this)) {}

void MeshComponent::WriteInternal(StructureSink* sink) const {
  if (texture_atlas_.Get() == nullptr) {
    sink->WritePod<int8>(false);
  } else {
    sink->WritePod<int8>(true);
    sink->WriteImage(texture_atlas_);
  }

  sink->WritePod<int64>(attribute_buf_.size());
  for (const VertexAttributes& attr : attribute_buf_) {
    sink->WritePoint(attr.position);
    sink->WritePoint(attr.tex_coord);
  }

  sink->WritePod<int64>(index_buf_.size());
  sink->WritePodArray(index_buf_.data(), index_buf_.size());

  sink->WritePoint(aabb_.GetMinPoint());
  sink->WritePoint(aabb_.GetMaxPoint());
}

REGISTER_SEURAT_COMPONENT(MeshComponent);

}  // namespace mesh
}  // namespace seurat
