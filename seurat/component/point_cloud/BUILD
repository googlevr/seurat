# Copyright 2017 Google Inc. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS-IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

load(
    "//seurat/dev:build_defs.bzl",
    "compiled_zipassets",
)

package_group(
    name = "friends",
    packages = [
        "//seurat/...",
    ],
)

package(
    default_visibility = [":friends"],
)

# Point cloud shader assets.
compiled_zipassets(
    name = "point_cloud_assets",
    srcs = [
        "point_cloud.frag",
        "point_cloud.vert",
    ],
    asset_file = "point_cloud_assets.iad",
)

cc_library(
    name = "point_cloud",
    srcs = [
        "point_cloud_component.cc",
        "point_cloud_renderable.cc",
    ],
    hdrs = [
        "point_cloud_component.h",
        "point_cloud_renderable.h",
    ],
    visibility = [":friends"],
    deps = [
        ":point_cloud_assets",
        "//seurat/base",
        "//seurat/component",
        "@com_google_absl//absl/types:span",
        "@ion//:graphicsmanager",
        "@ion//:ionbase",
        "@ion//:iongfx",
        "@ion//:iongfxutils",
        "@ion//:ionmath",
        "@ion//:ionremote",
        "@ion//:statetable",
    ],
    alwayslink = 1,
)
