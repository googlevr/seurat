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

#include "seurat/component/point_cloud/point_cloud_component.h"

#include <array>

#include "ion/math/matrix.h"
#include "ion/math/matrixutils.h"
#include "ion/math/transformutils.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/component/point_cloud/point_cloud_renderable.h"

namespace seurat {
namespace component {

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
    const std::vector<PointCloudComponent::VertexAttributes>& verts) {
  Range3f aabb;
  for (const auto& vert : verts) {
    aabb.ExtendByPoint(vert.position);
  }
  return aabb;
}

}  // namespace

// static
std::unique_ptr<PointCloudComponent> PointCloudComponent::Create(
    std::string label, std::vector<VertexAttributes> attribute_buf) {
  const auto aabb = ComputeBoundingBox(attribute_buf);
  return std::unique_ptr<PointCloudComponent>(new PointCloudComponent(
      std::move(label), std::move(attribute_buf), aabb));
}

// static
std::unique_ptr<const Component> PointCloudComponent::Create(
    std::string label, StructureSource* source) {
  std::vector<VertexAttributes> attribute_buf(source->ReadPod<int64>());
  for (VertexAttributes& attr : attribute_buf) {
    attr.position = source->ReadPoint<Point3f>();
  }

  const auto aabb_min = source->ReadPoint<Point3f>();
  const auto aabb_max = source->ReadPoint<Point3f>();
  const Range3f aabb(aabb_min, aabb_max);

  return std::unique_ptr<PointCloudComponent>(new PointCloudComponent(
      std::move(label), std::move(attribute_buf), aabb));
}

Renderable* PointCloudComponent::GetRenderable() const {
  return renderable_.get();
}

bool PointCloudComponent::operator==(const Component& other) const {
  if (!IsComponentEqualTo(other)) {
    return false;
  }
  const PointCloudComponent* other_mesh =
      static_cast<const PointCloudComponent*>(&other);
  if (other_mesh->GetAttributeBuffer() != GetAttributeBuffer()) {
    return false;
  }
  if (other_mesh->GetBoundingBox() != GetBoundingBox()) {
    return false;
  }
  return true;
}

PointCloudComponent::PointCloudComponent(
    std::string label, std::vector<VertexAttributes> attribute_buf,
    const ion::math::Range3f& aabb)
    : Component(std::move(label)),
      attribute_buf_(std::move(attribute_buf)),
      aabb_(aabb),
      renderable_(new PointCloudRenderable(this)) {}

void PointCloudComponent::WriteInternal(StructureSink* sink) const {
  sink->WritePod<int64>(attribute_buf_.size());
  for (const VertexAttributes& attr : attribute_buf_) {
    sink->WritePoint(attr.position);
  }

  sink->WritePoint(aabb_.GetMinPoint());
  sink->WritePoint(aabb_.GetMaxPoint());
}

REGISTER_SEURAT_COMPONENT(PointCloudComponent);

}  // namespace component
}  // namespace seurat
