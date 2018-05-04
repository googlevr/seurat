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
    "select_for_ion_platform",
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

base_sys_srcs = [
    "file_system_portable.cc",
]

base_sys_deps = []

cc_library(
    name = "base",
    srcs = [
        "bytestream.cc",
        "camera_util.cc",
        "file_system_common.cc",
        "ion_util_no_gl.cc",
        "math_util.cc",
        "parallel.cc",
        "progress.cc",
        "reporting.cc",
        "status.cc",
        "status_util.cc",
        "projective_camera.cc",
        "projective_camera_util.cc",
        "structured_io.cc",
    ] + base_sys_srcs,
    hdrs = [
        "array2d.h",
        "array2d_util.h",
        "array2d_view.h",
        "braced_init_list.h",
        "bytestream.h",
        "camera.h",
        "camera_util.h",
        "color.h",
        "file_system.h",
        "ion_util_no_gl.h",
        "math_util.h",
        "parallel.h",
        "progress.h",
        "projective_camera.h",
        "projective_camera_util.h",
        "reporting.h",
        "status.h",
        "status_util.h",
        "structured_io.h",
        "structured_io_tests_utils.h",
        "util.h",
    ],
    deps = [
        "//seurat/api:camera_portable_proto",
        "@ion//:ionbase",
        "@ion//:iongfx",
        "@ion//:iongfxutils",
        "@ion//:ionmath",
        "@ion//:ionport",
        "@ion//:ionprofile",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/base:core_headers",
        "@com_google_absl//absl/synchronization",
        "@ion//:ion_jsoncpp",
        "//seurat/api:math_portable_proto",
        "@com_google_absl//absl/types:span",
    ] + base_sys_deps,
)

cc_test(
    name = "projective_camera_test",
    size = "small",
    srcs = ["projective_camera_test.cc"],
    deps = [
        ":base",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "math_util_test",
    size = "small",
    srcs = ["math_util_test.cc"],
    deps = [
        ":base",
        "//seurat/testing",
        "//seurat/testing:test_flags",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "projective_camera_util_test",
    size = "small",
    srcs = ["projective_camera_util_test.cc"],
    deps = [
        ":base",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "camera_util_test",
    size = "small",
    srcs = ["camera_util_test.cc"],
    deps = [
        ":base",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "color_test",
    size = "small",
    srcs = ["color_test.cc"],
    deps = [
        ":base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "array2d_test",
    size = "small",
    srcs = ["array2d_test.cc"],
    deps = [
        ":base",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "array2d_util_test",
    size = "small",
    srcs = ["array2d_util_test.cc"],
    deps = [
        ":base",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "array2d_view_test",
    size = "small",
    srcs = ["array2d_view_test.cc"],
    deps = [
        ":base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "file_system_test",
    size = "small",
    srcs = [
        "file_system_test.cc",
    ],
    tags = [
        "no-android",
    ],
    deps = [
        ":base",
        "//seurat/testing",
        "//seurat/testing:test_flags",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionport",
    ],
)

cc_test(
    name = "ion_util_no_gl_test",
    size = "small",
    srcs = ["ion_util_no_gl_test.cc"],
    deps = [
        ":base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "parallel_test",
    size = "small",
    srcs = ["parallel_test.cc"],
    deps = [
        ":base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "util_test",
    size = "small",
    srcs = ["util_test.cc"],
    deps = [
        ":base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "reporting_test",
    size = "small",
    srcs = ["reporting_test.cc"],
    deps = [
        ":base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "status_test",
    size = "small",
    srcs = ["status_test.cc"],
    deps = [
        ":base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "braced_init_list_test",
    size = "small",
    srcs = ["braced_init_list_test.cc"],
    deps = [
        ":base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "status_util_test",
    size = "small",
    srcs = ["status_util_test.cc"],
    deps = [
        ":base",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "structured_io_test",
    size = "small",
    srcs = ["structured_io_test.cc"],
    deps = [
        ":base",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:iongfx",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "progress_test",
    size = "small",
    srcs = ["progress_test.cc"],
    deps = [
        ":base",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

test_suite(
    name = "base_tests",
    tests = [
        ":array2d_test",
        ":array2d_util_test",
        ":array2d_view_test",
        ":braced_init_list_test",
        ":camera_util_test",
        ":color_test",
        ":file_system_test",
        ":ion_util_no_gl_test",
        ":math_util_test",
        ":parallel_test",
        ":progress_test",
        ":projective_camera_test",
        ":projective_camera_util_test",
        ":reporting_test",
        ":status_test",
        ":structured_io_test",
        ":util_test",
    ],
)
