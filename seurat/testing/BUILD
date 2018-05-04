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
    name = "testing",
    testonly = 1,
    hdrs = [
        "ion_test_utils.h",
    ],
    deps = [
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_library(
    name = "test_flags",
    testonly = 1,
    srcs = [
        "test_flags.cc",
    ],
    hdrs = [
        "test_flags.h",
    ],
    deps = [
        "@ion//:ionbase",
    ],
)

cc_test(
    name = "ion_test_utils_test",
    size = "small",
    srcs = ["ion_test_utils_test.cc"],
    deps = [
        ":testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionbase",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "test_flags_test",
    size = "small",
    srcs = ["test_flags_test.cc"],
    deps = [
        ":test_flags",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionbase",
    ],
)

test_suite(
    name = "testing_tests",
    tests = [
        ":ion_test_utils_test",
        ":test_flags_test",
    ],
)
