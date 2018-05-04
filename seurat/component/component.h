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

#ifndef VR_SEURAT_COMPONENT_COMPONENT_H_
#define VR_SEURAT_COMPONENT_COMPONENT_H_

#include <memory>
#include <string>

#include "seurat/base/structured_io.h"
#include "seurat/component/renderable.h"

namespace seurat {
namespace component {

// An interface for components of a Seurat scene.
//
// Subclasses must be declared in their class declarations with
// DECLARE_SEURAT_COMPONENT(ClassName), and registered with
// REGISTER_SEURAT_COMPONENT(ClassName).
//
// Components represent a hierarchical scene resulting from processing raw data
// (e.g. raw layered depth images) into efficiently-renderable formats (e.g.
// planar tile geometry).
// For example, a scene may consist of a root component storing a collection of
// cube-map components, each containing up to 6 child components representing
// the faces, each of which may then have their own child components containing
// geometry generated from a layered depth image for that cube face.
//
// Components in turn generate Renderable instances, which represent the scene
// in a form suitable for rendering.
class Component {
 public:
  // Creates a Component by reading data from a StructureSource.  Subclasses
  // should provide their own static implementations of Create() with the same
  // signature.
  static std::unique_ptr<const Component> Create(base::StructureSource* source);

  // Writes a Component into a StructureSink.  This function is non-virtual; it
  // calls the WriteInternal() virtual function below.
  void Write(base::StructureSink* sink) const;

  // Allows implementations to initialize the label, or not.
  explicit Component(std::string label = {}) : label_(std::move(label)) {}
  virtual ~Component() = default;

  // Returns the unique type-id identifying the concrete subclass.  Sub-classes
  // should implement this function using the DECLARE_SEURAT_COMPONENT()
  // macro.  Note that as implemented by DECLARE_SEURAT_COMPONENT(), this
  // function returns a string literal, which is defined to have static storage
  // duration; hence it is safe to pointer-compare results of GetTypeId()
  // directly.
  virtual const char* GetTypeId() const = 0;

  // Comparison operator.  Useful for testing.
  virtual bool operator==(const Component& other) const = 0;
  bool operator!=(const Component& other) const { return !(*this == other); }

  // Returns a pointer to this component's Renderable.
  virtual Renderable* GetRenderable() const = 0;

  // Returns the label.
  const std::string& GetLabel() const { return label_; }

  // Changes the component label to |label|.
  void SetLabel(const std::string& label) { label_ = label; }

 protected:
  // Writes a Component's internal representation into a base::StructureSink.
  //
  // NOTE! Component implementations are responsible for serializing the label.
  virtual void WriteInternal(base::StructureSink* sink) const = 0;

  // Tests for Component type and member equality; assists in making equality
  // operator implementations, while leaving operator == a pure method. Use this
  // in implementations.
  bool IsComponentEqualTo(const Component& other) const {
    return GetTypeId() == other.GetTypeId() && label_ == other.label_;
  }

 private:
  // Holds a human-readable name or label for the Component, as well as any
  // Renderable data it generates. Labels should follow google variable naming
  // conventions and are accessed programmatically in some cases.
  std::string label_;

  template <typename T>
  friend class ComponentRegisterer;

  // Registers a factory function for creating a Component.  For internal use;
  // use only through REGISTER_SEURAT_COMPONENT().  Not thread safe!
  static void RegisterFactory(
      const char* name,
      std::unique_ptr<const Component> (*func)(std::string label,
                                               base::StructureSource* source));
};

// Declare a Component.  Use DECLARE_SEURAT_COMPONENT in the Component
// subclass's class declaration.  This implementation of GetTypeId() returns a
// string literal, which is defined to have static storage duration; hence it is
// safe to compare results of GetTypeId() directly.
#define DECLARE_SEURAT_COMPONENT(ClassName)                \
  friend class component::ComponentRegisterer<ClassName>;     \
  static const char* GetStaticTypeId() { return #ClassName; } \
  const char* GetTypeId() const override { return GetStaticTypeId(); }

// Registers a Component.  Use REGISTER_SEURAT_COMPONENT once in a
// translation unit containing the Component subclass's implementation.
#define REGISTER_SEURAT_COMPONENT(ClassName)  \
  namespace {                                    \
  component::ComponentRegisterer<ClassName>      \
      s_##ClassName##_creation_factory_instance; \
  }
// Templated registration class for use with REGISTER_SEURAT_COMPONENT().
template <typename T>
class ComponentRegisterer {
 public:
  ComponentRegisterer() {
    component::Component::RegisterFactory(T::GetStaticTypeId(), &T::Create);
  }
};

}  // namespace component
}  // namespace seurat

#endif  // VR_SEURAT_COMPONENT_COMPONENT_H_
