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

#ifndef VR_SEURAT_COMPONENT_POINT_CLOUD_POINT_CLOUD_COMPONENT_H_
#define VR_SEURAT_COMPONENT_POINT_CLOUD_POINT_CLOUD_COMPONENT_H_

#include <array>

#include "ion/gfx/image.h"
#include "ion/math/vector.h"
#include "seurat/component/component.h"

namespace seurat {
namespace component {

// A component for rendering a point cloud.
class PointCloudComponent : public component::Component {
 public:
  struct VertexAttributes {
    VertexAttributes() = default;
    explicit VertexAttributes(const ion::math::Point3f& position)
        : position(position) {}

    bool operator==(const VertexAttributes& rhs) const {
      return position == rhs.position;
    }

    ion::math::Point3f position;
  };

  ~PointCloudComponent() override = default;

  // Creates a new PointCloudComponent with the given attribute buffer.
  static std::unique_ptr<PointCloudComponent> Create(
      std::string label, std::vector<VertexAttributes> attribute_buf);

  // Component implementation.
  static std::unique_ptr<const component::Component> Create(
      std::string label, base::StructureSource* source);
  component::Renderable* GetRenderable() const override;
  bool operator==(const component::Component& other) const override;

  // Returns true if this point cloud component contains no points.
  bool IsEmpty() const { return attribute_buf_.empty(); }

  // Returns the attribute buffer for this point cloud, suitable for sending
  // directly to OpenGL.
  const std::vector<VertexAttributes>& GetAttributeBuffer() const {
    return attribute_buf_;
  }

  // Returns an axis-aligned bounding box for the points.
  const ion::math::Range3f& GetBoundingBox() const { return aabb_; }

 private:
  PointCloudComponent(std::string label,
                      std::vector<VertexAttributes> attribute_buf,
                      const ion::math::Range3f& aabb);

  // Component implementation.
  void WriteInternal(base::StructureSink* sink) const override;

  // Interleaved attribute-buffer which can be sent directly to opengl for
  // rendering.
  const std::vector<VertexAttributes> attribute_buf_;

  // Axis-aligned bounding box for the points.
  const ion::math::Range3f aabb_;

  // The renderable used for rendering this component.
  const std::unique_ptr<component::Renderable> renderable_;

  DECLARE_SEURAT_COMPONENT(PointCloudComponent);
};

}  // namespace component
}  // namespace seurat

#endif  // VR_SEURAT_COMPONENT_POINT_CLOUD_POINT_CLOUD_COMPONENT_H_
