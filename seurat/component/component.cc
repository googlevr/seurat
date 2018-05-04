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

#include "seurat/component/component.h"

#include <unordered_map>

#include "ion/base/staticsafedeclare.h"

namespace seurat {
namespace component {
namespace {

using FactoryFunc = std::unique_ptr<const Component> (*)(
    std::string label, base::StructureSource* source);
using FactoryMap = std::unordered_map<std::string, FactoryFunc>;

FactoryMap* GetFactoryMap() {
  ION_DECLARE_SAFE_STATIC_POINTER(FactoryMap, s_factory_map);
  return s_factory_map;
}

}  // namespace

// static
std::unique_ptr<const Component> Component::Create(
    base::StructureSource* source) {
  const std::string name = source->ReadString();
  std::string label = source->ReadString();
  const FactoryMap* const factory_map = GetFactoryMap();
  const auto iter = factory_map->find(name);
  if (iter == factory_map->end()) {
    LOG(ERROR) << "Component name=\"" << name << "\" not found for creation";
    return nullptr;
  }
  return iter->second(std::move(label), source);
}

void Component::Write(base::StructureSink* sink) const {
  sink->WriteString(GetTypeId());
  sink->WriteString(label_);
  WriteInternal(sink);
}

// static
void Component::RegisterFactory(const char* name, FactoryFunc func) {
  FactoryMap* const factory_map = GetFactoryMap();
  factory_map->insert(std::make_pair(std::string(name), func));
}

}  // namespace component
}  // namespace seurat
