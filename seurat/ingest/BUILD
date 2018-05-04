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

cc_library(
    name = "ingest",
    srcs = [
        "clamping_view_group_loader.cc",
        "json_converter.cc",
        "json_utils.cc",
        "json_validator.cc",
        "ldi_loader.cc",
        "pattern_loader.cc",
        "point_cloud_assembler.cc",
        "point_cloud_assembler_factory.cc",
        "proto_view_group_loader.cc",
        "single_face_view_group_loader.cc",
        "sphere_view_group_loader.cc",
        "subset_view_group_loader.cc",
        "tone_mapping_view_group_loader.cc",
        "view_group_loader_factory.cc",
        "view_group_loader_util.cc",
    ],
    hdrs = [
        "clamping_view_group_loader.h",
        "json_converter.h",
        "json_utils.h",
        "json_validator.h",
        "ldi_loader.h",
        "pattern_loader.h",
        "point_cloud_assembler.h",
        "point_cloud_assembler_factory.h",
        "proto_view_group_loader.h",
        "single_face_view_group_loader.h",
        "sphere_view_group_loader.h",
        "subset_view_group_loader.h",
        "tone_mapping_view_group_loader.h",
        "view_group_loader.h",
        "view_group_loader_factory.h",
        "view_group_loader_util.h",
    ],
    deps = [
        "//seurat/api:api_portable_proto",
        "//seurat/api:camera_portable_proto",
        "//seurat/api:image_portable_proto",
        "//seurat/api:math_portable_proto",
        "//seurat/base",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/image:exr_io",
        "@com_google_absl//absl/strings",
        "@ion//:ion_jsoncpp",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

cc_library(
    name = "ingest_test_utils",
    testonly = 1,
    srcs = [
        "ldi_loader_test_utils.cc",
        "view_group_loader_test_utils.cc",
    ],
    hdrs = [
        "ldi_loader_test_utils.h",
        "view_group_loader_test_utils.h",
    ],
    data = [
        ":testdata/arnold_deep_exr.exr",
    ],
    deps = [
        ":ingest",
        "//seurat/base",
        "//seurat/geometry",
        "//seurat/image",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "view_group_loader_factory_test",
    size = "small",
    srcs = ["view_group_loader_factory_test.cc"],
    data = [
        ":testdata/color.exr",
        ":testdata/depth.exr",
        ":testdata/manifest_v3.json",
    ],
    deps = [
        ":ingest",
        "//seurat/base",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/testing",
        "//seurat/testing:test_flags",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "view_group_loader_util_test",
    size = "small",
    srcs = ["view_group_loader_util_test.cc"],
    deps = [
        ":ingest",
        ":ingest_test_utils",
        "//seurat/base",
        "//seurat/image",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "json_utils_test",
    size = "small",
    srcs = ["json_utils_test.cc"],
    data = [
        ":testdata/invalid_manifest.json",
        ":testdata/manifest_v3.json",
    ],
    deps = [
        ":ingest",
        ":ingest_test_utils",
        "//seurat/base",
        "//seurat/testing",
        "//seurat/testing:test_flags",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ion_jsoncpp",
    ],
)

cc_test(
    name = "json_validator_test",
    size = "small",
    srcs = ["json_validator_test.cc"],
    data = [
        ":testdata/manifest_v3.json",
    ],
    deps = [
        ":ingest",
        ":ingest_test_utils",
        "//seurat/base",
        "//seurat/testing",
        "//seurat/testing:test_flags",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ion_jsoncpp",
    ],
)

cc_test(
    name = "json_converter_test",
    size = "small",
    srcs = ["json_converter_test.cc"],
    data = [
        ":testdata/manifest_v3.json",
    ],
    deps = [
        ":ingest",
        ":ingest_test_utils",
        "//seurat/api:api_portable_proto",
        "//seurat/api:camera_portable_proto",
        "//seurat/api:image_portable_proto",
        "//seurat/api:math_portable_proto",
        "//seurat/base",
        "//seurat/testing",
        "//seurat/testing:test_flags",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ion_jsoncpp",
    ],
)

cc_test(
    name = "ldi_loader_test",
    size = "small",
    srcs = ["ldi_loader_test.cc"],
    data = [
        ":testdata/arnold_deep_exr.exr",
        ":testdata/depth_image_with_transparency.exr",
    ],
    deps = [
        ":ingest",
        ":ingest_test_utils",
        "//seurat/api:image_portable_proto",
        "//seurat/base",
        "//seurat/image",
        "//seurat/testing",
        "//seurat/testing:test_flags",
        "@com_google_absl//absl/memory",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "proto_view_group_loader_test",
    size = "small",
    srcs = ["proto_view_group_loader_test.cc"],
    deps = [
        ":ingest",
        ":ingest_test_utils",
        "//seurat/api:api_portable_proto",
        "//seurat/api:math_portable_proto",
        "//seurat/base",
        "//seurat/image",
        "//seurat/image:image_test_utils",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

cc_test(
    name = "sphere_view_group_loader_test",
    size = "small",
    srcs = ["sphere_view_group_loader_test.cc"],
    deps = [
        ":ingest",
        ":ingest_test_utils",
        "//seurat/base",
        "//seurat/image",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

cc_test(
    name = "subset_view_group_loader_test",
    size = "small",
    srcs = ["subset_view_group_loader_test.cc"],
    deps = [
        ":ingest",
        ":ingest_test_utils",
        "//seurat/base",
        "//seurat/image",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

cc_test(
    name = "clamping_view_group_loader_test",
    size = "small",
    srcs = ["clamping_view_group_loader_test.cc"],
    deps = [
        ":ingest",
        ":ingest_test_utils",
        "//seurat/base",
        "//seurat/image",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

cc_test(
    name = "single_face_view_group_loader_test",
    size = "small",
    srcs = ["single_face_view_group_loader_test.cc"],
    deps = [
        ":ingest",
        ":ingest_test_utils",
        "//seurat/base",
        "//seurat/image",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

cc_test(
    name = "pattern_loader_test",
    size = "small",
    srcs = ["pattern_loader_test.cc"],
    deps = [
        ":ingest",
        ":ingest_test_utils",
        "//seurat/base",
        "//seurat/image",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

cc_test(
    name = "point_cloud_assembler_test",
    size = "small",
    srcs = ["point_cloud_assembler_test.cc"],
    deps = [
        ":ingest",
        ":ingest_test_utils",
        "//seurat/base",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

test_suite(
    name = "ingest_tests",
    tests = [
        "view_group_loader_util_test",
        ":clamping_view_group_loader_test",
        ":json_converter_test",
        ":json_utils_test",
        ":json_validator_test",
        ":ldi_loader_test",
        ":pattern_loader_test",
        ":point_cloud_assembler_test",
        ":proto_view_group_loader_test",
        ":single_face_view_group_loader_test",
        ":sphere_view_group_loader_test",
        ":subset_view_group_loader_test",
        ":view_group_loader_factory_test",
    ],
)
