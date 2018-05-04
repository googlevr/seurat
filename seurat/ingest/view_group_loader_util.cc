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

#include "seurat/ingest/view_group_loader_util.h"

#include <future>  // NOLINT(build/c++11)
#include <list>

#include "ion/math/transformutils.h"
#include "seurat/base/reporting.h"

namespace seurat {
namespace ingest {

ion::math::Range3f ComputeHeadbox(const ViewGroupLoader& view_group_loader) {
  ion::math::Range3f headbox;
  std::vector<std::shared_ptr<base::Camera>> cameras;
  for (int i = 0; i < view_group_loader.GetNumViewGroups(); ++i) {
    // TODO(ernstm): Propagate this status up.
    base::Status status = view_group_loader.LoadViewGroup(i, &cameras, nullptr);
    if (!status.ok()) {
      base::SeuratFatal(status.error_message());
    }
    for (const auto& camera : cameras) {
      const ion::math::Point3f position_world =
          camera->GetWorldFromEye() * ion::math::Point3f::Zero();
      headbox.ExtendByPoint(position_world);
    }
  }
  return headbox;
}

namespace {

// The result of loading a view group.
struct ViewGroupResult {
  std::vector<std::shared_ptr<base::Camera>> cameras;
  std::vector<image::Ldi4f> ldis;
  base::Status status;
};

std::future<ViewGroupResult> LoadAsync(const ViewGroupLoader& loader,
                                       int index) {
  // Note that this grabs a reference to |loader|, so the caller must be careful
  // to wait until all threads have completed before moving on.
  return std::async(std::launch::async, [&loader, index]() {
    ViewGroupResult group;
    group.status = loader.LoadViewGroup(index, &group.cameras, &group.ldis);
    return group;
  });
}

}  // namespace

base::Status ForEachViewGroupPrefetching(
    const ViewGroupLoader& loader,
    const std::function<base::Status(std::vector<std::shared_ptr<base::Camera>>,
                                     std::vector<image::Ldi4f>)>& func) {
  const int kNumInFlightLoads = 4;

  int num_groups = loader.GetNumViewGroups();

  // Pending view groups which have started/finished loading and are to be
  // processed.
  std::list<std::future<ViewGroupResult>> groups;

  // Start loading the first couple of groups.
  //
  // Loading multiple here ensures that there is some loading which can overlap
  // with the processing performed by |func|.
  int next_group_to_fetch = 0;
  for (int i = 0; i < kNumInFlightLoads; ++i) {
    if (next_group_to_fetch < num_groups) {
      groups.push_back(LoadAsync(loader, next_group_to_fetch));
      next_group_to_fetch++;
    }
  }

  // Loop over all groups, processing them as they become available & kicking
  // off additional async loads down the line.
  for (int i = 0; i < num_groups; ++i) {
    groups.front().wait();
    ViewGroupResult result = groups.front().get();
    groups.pop_front();
    if (!result.status.ok()) {
      for (const auto& future : groups) {
        // Wait for all (other) outstanding async operations to complete before
        // returning, since those threads may be holding a reference to
        // |loader|.
        future.wait();
      }
      return result.status;
    }

    if (next_group_to_fetch < num_groups) {
      // Start prefetching another group.
      groups.push_back(LoadAsync(loader, next_group_to_fetch));
      next_group_to_fetch++;
    }

    SEURAT_RETURN_IF_ERROR(
        func(std::move(result.cameras), std::move(result.ldis)));
  }
  return base::OkStatus();
}

}  // namespace ingest
}  // namespace seurat
