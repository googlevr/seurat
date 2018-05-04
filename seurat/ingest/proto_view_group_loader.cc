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

#include "seurat/ingest/proto_view_group_loader.h"

#include <memory>

#include "ion/math/transformutils.h"
#include "absl/strings/str_cat.h"
#include "seurat/base/camera_util.h"
#include "seurat/base/math_util.h"
#include "seurat/base/parallel.h"
#include "seurat/base/status_util.h"
#include "seurat/base/util.h"

namespace seurat {
namespace ingest {
namespace {

// Returns the filename(s) referenced in the |ldi_proto| in a string. This is
// useful to compose error messages.
std::string GetLdiFilenames(const api::proto::Ldi& ldi_proto) {
  switch (ldi_proto.ldi_case()) {
    case api::proto::Ldi::LdiCase::LDI_NOT_SET:
      return "unspecified file";
    case api::proto::Ldi::LdiCase::kLdiFile:
      return ldi_proto.ldi_file().path();
    case api::proto::Ldi::LdiCase::kDepthImageFile:
      return ldi_proto.depth_image_file().color().path() + " and " +
             ldi_proto.depth_image_file().color().path();
  }
}

}  // namespace

ProtoViewGroupLoader::ProtoViewGroupLoader(
    const api::proto::Capture& capture, int thread_count,
    std::unique_ptr<LdiLoader> ldi_loader,
    std::shared_ptr<base::FileSystem> file_system)
    : capture_(capture),
      thread_count_(thread_count),
      ldi_loader_(std::move(ldi_loader)),
      file_system_(std::move(file_system)) {}

int ProtoViewGroupLoader::GetNumViewGroups() const {
  return capture_.view_groups_size();
}

base::Status ProtoViewGroupLoader::LoadViewGroup(
    int view_group_index, std::vector<std::shared_ptr<base::Camera>>* cameras,
    std::vector<image::Ldi4f>* ldis) const {
  if (view_group_index >= GetNumViewGroups()) {
    return base::OutOfRangeError(absl::StrCat(
        "view_group_index ", view_group_index,
        " out of range. Number of view groups is ", GetNumViewGroups()));
  }

  api::proto::ViewGroup view_group = capture_.view_groups(view_group_index);
  if (ldis != nullptr) {
    ldis->clear();
    ldis->resize(view_group.views_size());
  }
  if (cameras != nullptr) {
    cameras->clear();
    cameras->resize(view_group.views_size());
  }
  std::vector<base::Status> per_view_status(view_group.views_size());

  // Note, protobuf supports multiple reader threads, provided there are no
  // writers.
  base::ParallelFor(
      thread_count_, view_group.views_size(), [&](int view_index) {
        api::proto::View view = view_group.views(view_index);

        std::unique_ptr<base::Camera> camera;
        if (capture_.has_headbox_center()) {
          const ion::math::Vector3f translation =
              ion::math::Point3f::Zero() -
              base::Point3fFromProto(capture_.headbox_center());
          camera = base::TranslateCamera(translation,
                                         base::CameraFromProto(view.camera()));
        } else {
          camera = base::CameraFromProto(view.camera());
        }

        ion::math::Vector2i image_size = camera->GetImageSize();

        if (cameras != nullptr) {
          (*cameras)[view_index] = std::move(camera);
        }

        if (ldis != nullptr) {
          api::proto::Ldi ldi_proto = view.ldi();

          // Try loading the ldi. Set status and early out if it fails.
          image::Ldi4f ldi;
          per_view_status[view_index] =
              ldi_loader_->Load(ldi_proto, file_system_.get(), &ldi);
          if (!per_view_status[view_index].ok()) {
            return;
          }

          // Verify that the size of the LDI matches the image size specified in
          // the proto. Set status and early out if they don't match.
          ion::math::Vector2i ldi_size = ldi.GetSize();
          if (ldi_size != image_size) {
            per_view_status[view_index] = base::InvalidArgumentError(
                absl::StrCat("Actual image size ", ldi_size[0], "x",
                             ldi_size[1], " of ", GetLdiFilenames(ldi_proto),
                             " does not match image size ", image_size[0], "x",
                             image_size[1], " specified in manifest file."));
            return;
          }

          // Store the LDI only if there were no errors.
          if (per_view_status[view_index].ok()) {
            ldis->at(view_index) = std::move(ldi);
          }
        }
      });

  // If an error occurred, return only the status from the first thread that
  // failed.
  for (auto const& view_status : per_view_status) {
    SEURAT_RETURN_IF_ERROR(view_status);
  }

  return base::OkStatus();
}

}  // namespace ingest
}  // namespace seurat
