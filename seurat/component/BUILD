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

# Shape shader assets.
compiled_zipassets(
    name = "shape_assets",
    srcs = glob([
        "shaders/shape.*",
    ]),
    asset_file = "shape_assets.iad",
)

cc_library(
    name = "component",
    srcs = [
        "component.cc",
        "group_component.cc",
        "ice_io.cc",
        "renderable_util.cc",
        "shape_component.cc",
        "shape_component_utils.cc",
        "shape_renderable.cc",
        "transform_component.cc",
    ],
    hdrs = [
        "component.h",
        "group_component.h",
        "ice_io.h",
        "renderable.h",
        "renderable_util.h",
        "shape_component.h",
        "shape_component_utils.h",
        "shape_renderable.h",
        "transform_component.h",
    ],
    deps = [
        ":shape_assets",
        "//seurat/base",
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

cc_library(
    name = "component_tests_utils",
    testonly = 1,
    srcs = ["component_tests_utils.cc"],
    hdrs = ["component_tests_utils.h"],
    deps = [
        ":component",
        "//seurat/base",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "component_io_test",
    size = "small",
    srcs = ["component_io_test.cc"],
    deps = [
        ":component",
        ":component_tests_utils",
        "//seurat/base",
        "//seurat/image",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:iongfx",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "ice_io_test",
    size = "small",
    srcs = ["ice_io_test.cc"],
    deps = [
        ":component",
        ":component_tests_utils",
        "//seurat/base",
        "//seurat/image",
        "//seurat/testing",
        "//seurat/testing:test_flags",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:iongfx",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "component_test",
    size = "small",
    srcs = ["component_test.cc"],
    deps = [
        ":component",
        ":component_tests_utils",
        "//seurat/base",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "group_component_test",
    size = "small",
    srcs = ["group_component_test.cc"],
    deps = [
        ":component",
        ":component_tests_utils",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "transform_component_test",
    size = "small",
    srcs = ["transform_component_test.cc"],
    deps = [
        ":component",
        ":component_tests_utils",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "shape_component_test",
    size = "small",
    srcs = ["shape_component_test.cc"],
    deps = [
        ":component",
        ":component_tests_utils",
        "//seurat/base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "shape_component_utils_test",
    size = "small",
    srcs = ["shape_component_utils_test.cc"],
    deps = [
        ":component",
        ":component_tests_utils",
        "//seurat/base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

test_suite(
    name = "component_tests",
    tests = [
        ":component_io_test",
        ":component_test",
        ":group_component_test",
        ":ice_io_test",
        ":shape_component_test",
        ":shape_component_utils_test",
        ":transform_component_test",
    ],
)
