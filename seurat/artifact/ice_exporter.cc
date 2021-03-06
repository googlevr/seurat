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

#include "seurat/artifact/ice_exporter.h"

#include "absl/strings/str_cat.h"
#include "seurat/base/status.h"
#include "seurat/base/structured_io.h"
#include "seurat/component/ice_io.h"

namespace seurat {
namespace artifact {

using base::StructureSink;
using component::Component;

// Exporter implementation.
base::Status IceExporter::Process(Artifact* artifact) const {
  if (!artifact->component) {
    return base::UnimplementedError("No component to export.");
  }
  std::string filename = absl::StrCat(basename_, ".ice");
  std::string encoded;
  base::StringByteSink sink(&encoded);
  StructureSink component_sink(&sink);
  component::WriteIce(*artifact->component, &component_sink);

  return filesystem_->SetContents(filename, encoded);
}

}  // namespace artifact
}  // namespace seurat
