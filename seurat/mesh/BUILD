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

# Mesh component shader assets.
compiled_zipassets(
    name = "mesh_component_assets",
    srcs = [
        "shaders/mesh.frag",
        "shaders/mesh.vert",
        "shaders/solid_depth_mesh.frag",
        "shaders/solid_depth_mesh.vert",
    ],
    asset_file = "mesh_component_assets.iad",
)

cc_library(
    name = "mesh",
    srcs = [
        "mesh_component.cc",
        "mesh_renderable.cc",
        "mesh_renderable_util.cc",
    ],
    hdrs = [
        "mesh_component.h",
        "mesh_renderable.h",
        "mesh_renderable_util.h",
    ],
    deps = [
        ":mesh_component_assets",
        "//seurat/base",
        "//seurat/component",
        "//seurat/geometry",
        "//seurat/image",
        "@com_google_absl//absl/types:span",
        "@ion//:graphicsmanager",
        "@ion//:ionbase",
        "@ion//:iongfx",
        "@ion//:iongfxutils",
        "@ion//:ionmath",
        "@ion//:statetable",
    ],
    alwayslink = 1,
)

cc_library(
    name = "mesh_component_util",
    srcs = [
        "mesh_component_util.cc",
    ],
    hdrs = [
        "mesh_component_util.h",
    ],
    deps = [
        ":mesh",
        "//seurat/base",
        "//seurat/component",
        "//seurat/geometry",
        "//seurat/image",
        "@ion//:iongfx",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "mesh_component_test",
    size = "small",
    srcs = ["mesh_component_test.cc"],
    deps = [
        ":mesh",
        "//seurat/base",
        "//seurat/component:component_tests_utils",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "mesh_component_util_test",
    size = "small",
    srcs = ["mesh_component_util_test.cc"],
    deps = [
        ":mesh_component_util",
        "//seurat/base",
        "//seurat/geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

test_suite(
    name = "mesh_tests",
    tests = [
        ":mesh_component_test",
        ":mesh_component_util_test",
    ],
)
