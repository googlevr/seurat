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
    name = "rgba",
    srcs = [
        "rgba_compressor_util.cc",
        "rgba_encoding_selector.cc",
        "rgba_nop_compressor.cc",
        "rgba_rate_resizer.cc",
        "rgba_selecting_compressor.cc",
    ],
    hdrs = [
        "rgba_compressor.h",
        "rgba_compressor_util.h",
        "rgba_encoding_selector.h",
        "rgba_nop_compressor.h",
        "rgba_rate_resizer.h",
        "rgba_selecting_compressor.h",
    ],
    deps = [
        "//seurat/base",
        "//seurat/compressor",
        "//seurat/compressor/image_metrics",
        "//seurat/compressor/resampler",
        "//seurat/image",
        "@com_google_absl//absl/types:span",
        "@ion//:ionmath",
        "@ion//:ionprofile",
    ],
)

cc_test(
    name = "rgba_rate_resizer_test",
    size = "small",
    srcs = ["rgba_rate_resizer_test.cc"],
    deps = [
        ":rgba",
        "//seurat/base",
        "//seurat/compressor/image_metrics",
        "//seurat/compressor/resampler",
        "//seurat/image",
        "//seurat/image:image_test_utils",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "rgba_encoding_selector_test",
    size = "small",
    srcs = ["rgba_encoding_selector_test.cc"],
    deps = [
        ":rgba",
        "//seurat/base",
        "//seurat/compressor/image_metrics",
        "//seurat/compressor/resampler",
        "//seurat/image",
        "//seurat/image:image_test_utils",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "rgba_selecting_compressor_test",
    size = "small",
    srcs = ["rgba_selecting_compressor_test.cc"],
    data = [
    ],
    deps = [
        ":rgba",
        "//seurat/base",
        "//seurat/compressor/image_metrics",
        "//seurat/compressor/resampler",
        "//seurat/image",
        "//seurat/image:image_test_utils",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

test_suite(
    name = "rgba_tests",
    tests = [
        ":rgba_encoding_selector_test",
        ":rgba_rate_resizer_test",
        ":rgba_selecting_compressor_test",
    ],
)
