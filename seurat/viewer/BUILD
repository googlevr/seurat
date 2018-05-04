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
    "init_ion_settings",
    "cc_test",
    "test_suite",
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

# Required for compiled_zipassets()
init_ion_settings()

components = [
    "//seurat/component/point_cloud",
    "//seurat/mesh",
]

cc_library(
    name = "viewer",
    srcs = [
        "rendering_pipeline.cc",
        "scene.cc",
        "viewer_util.cc",
    ],
    hdrs = [
        "rendering_pipeline.h",
        "scene.h",
        "viewer_util.h",
    ],
    deps = [
        "@ion//:ionbase",
        "@ion//:graphicsmanager",
        "@ion//:statetable",
        "@ion//:iongfxutils",
        "@ion//:ionimage",
        "@ion//:iongfx",
        "@ion//:ionmath",
        "@ion//:ionportgfx",
        "@com_google_absl//absl/strings",
        "//seurat/component",
    ] + components,
)

cc_test(
    name = "scene_test",
    size = "small",
    srcs = ["scene_test.cc"],
    deps = [
        ":viewer",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:iongfx",
        "@ion//:ionmath",
    ],
)

test_suite(
    name = "viewer_tests",
    tests = [
        ":scene_test",
    ],
)
