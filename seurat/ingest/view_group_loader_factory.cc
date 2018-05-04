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

#include "seurat/ingest/view_group_loader_factory.h"

#include <string>

#include "json/json.h"
#include "seurat/api/api.pb.h"
#include "seurat/base/projective_camera.h"
#include "seurat/base/projective_camera_util.h"
#include "seurat/base/status.h"
#include "seurat/ingest/clamping_view_group_loader.h"
#include "seurat/ingest/json_converter.h"
#include "seurat/ingest/json_utils.h"
#include "seurat/ingest/json_validator.h"
#include "seurat/ingest/ldi_loader.h"
#include "seurat/ingest/proto_view_group_loader.h"
#include "seurat/ingest/single_face_view_group_loader.h"
#include "seurat/ingest/subset_view_group_loader.h"

namespace seurat {
namespace ingest {

// static
base::Status ViewGroupLoaderFactory::CreateViewGroupLoader(
    absl::string_view manifest_filename,
    std::unique_ptr<ViewGroupLoader>* view_group_loader) const {
  if (!parameters_.pattern_parameters.name.empty()) {
    view_group_loader->reset(new PatternLoader(parameters_.pattern_parameters));
  }

  base::Status status;
  std::string manifest_string;
  SEURAT_RETURN_IF_ERROR(
      file_system_->GetContents(manifest_filename, &manifest_string));

  api::proto::Capture capture;

  // Check if this is a JSON manifest. If so, convert to Proto. Otherwise load
  // the proto from the file contents.
  Json::Value root;
  status = JsonUtils::ParseJsonFromString(manifest_string, &root);
  if (status.ok()) {
    SEURAT_RETURN_IF_ERROR(JsonValidator::ValidateCapture(root));
    JsonConverter::ConvertCapture(root, &capture);
  } else {
    if (!capture.ParseFromString(manifest_string)) {
      return base::InvalidArgumentError("Invalid manifest.");
    }
  }

  // A filesystem with a root path containing the manifest file.
  auto ldi_file_system = std::make_shared<base::FileSystem>(
      file_system_->SplitPath(file_system_->GetAbsolutePath(manifest_filename))
          .first);

  std::unique_ptr<LdiLoader> ldi_loader(new LdiLoader());
  view_group_loader->reset(
      new ProtoViewGroupLoader(capture, parameters_.thread_count,
                               std::move(ldi_loader), ldi_file_system));

  view_group_loader->reset(new ingest::SubsetViewGroupLoader(
      std::move(*view_group_loader), parameters_.max_views));

  if (parameters_.skybox_radius > 0.0) {
    // A depth of 0 indicates infinity.
    //
    // TODO(b/36357462):  Remove this once the ViewGroupLoader can read this
    // flag from the camera parameters in the json/proto.
    const bool kZeroIsInfinite = true;
    view_group_loader->reset(new ingest::ClampingViewGroupLoader(
        parameters_.skybox_radius, kZeroIsInfinite,
        std::move(*view_group_loader)));
  }
  if (!parameters_.single_face.empty()) {
    geometry::CubeFace face =
        geometry::CubeFaceFromString(parameters_.single_face);
    view_group_loader->reset(new ingest::SingleFaceViewGroupLoader(
        face, parameters_.thread_count, std::move(*view_group_loader)));
  }

  return base::OkStatus();
}

}  // namespace ingest
}  // namespace seurat
