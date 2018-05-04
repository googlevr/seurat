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
    "portable_proto_library",
    "cc_test",
    "test_suite",
)

package_group(
    name = "friends",
    packages = [
        "//seurat/...",
    ],
)

package(default_visibility = [":friends"])

proto_library(
    name = "api_proto",
    srcs = [
        "api.proto",
    ],
    deps = [
        ":camera_proto",
        ":image_proto",
        ":math_proto",
    ],
)

cc_proto_library(
    name = "api_portable_proto",
    deps = [":api_proto"],
)

proto_library(
    name = "math_proto",
    srcs = [
        "math.proto",
    ],
)

cc_proto_library(
    name = "math_portable_proto",
    deps = [":math_proto"],
)

proto_library(
    name = "camera_proto",
    srcs = [
        "camera.proto",
    ],
    deps = [":math_proto"],
)

cc_proto_library(
    name = "camera_portable_proto",
    deps = [":camera_proto"],
)

proto_library(
    name = "image_proto",
    srcs = [
        "image.proto",
    ],
)

cc_proto_library(
    name = "image_portable_proto",
    deps = [":image_proto"],
)

cc_library(
    name = "api",
    srcs = [
        "ldi_exporter.cc",
    ],
    hdrs = [
        "ldi_exporter.h",
    ],
    deps = [
        "//seurat/base",
        "//seurat/image",
        "//seurat/image:exr_io",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "ldi_exporter_test",
    size = "small",
    srcs = ["ldi_exporter_test.cc"],
    deps = [
        ":api",
        "//seurat/image:exr_io",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

test_suite(
    name = "api_tests",
    tests = [
        ":ldi_exporter_test",
    ],
)
