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
    name = "selection",
    srcs = [
        "bisection_solver.cc",
        "cost_calculator.cc",
        "parallel_solver.cc",
        "relaxed_solver.cc",
        "selection_problem.cc",
        "selection_util.cc",
        "simplifying_solver.cc",
        "solver_factory.cc",
        "subgradient_descent_solver.cc",
    ],
    hdrs = [
        "bisection_solver.h",
        "cost_calculator.h",
        "parallel_solver.h",
        "relaxed_solver.h",
        "selection_problem.h",
        "selection_solver.h",
        "selection_util.h",
        "sequential_solver.h",
        "simplifying_solver.h",
        "solver_factory.h",
        "subgradient_descent_solver.h",
    ],
    deps = [
        "//seurat/base",
        "@com_google_absl//absl/types:span",
        "@ion//:ionbase",
    ],
)

cc_library(
    name = "selection_test_util",
    testonly = 1,
    srcs = [
        "selection_test_util.cc",
    ],
    hdrs = [
        "selection_test_util.h",
    ],
    deps = [
        ":selection",
        "//seurat/base",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "bisection_solver_test",
    size = "small",
    srcs = [
        "bisection_solver_test.cc",
    ],
    tags = [
        "no-android",
    ],
    deps = [
        ":selection",
        ":selection_test_util",
        "//seurat/base",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "cost_calculator_test",
    size = "small",
    srcs = [
        "cost_calculator_test.cc",
    ],
    tags = [
        "no-android",
    ],
    deps = [
        ":selection",
        "//seurat/base",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "parallel_solver_test",
    size = "small",
    srcs = [
        "parallel_solver_test.cc",
    ],
    tags = [
        "no-android",
    ],
    deps = [
        ":selection",
        ":selection_test_util",
        "//seurat/base",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "relaxed_solver_test",
    size = "small",
    srcs = [
        "relaxed_solver_test.cc",
    ],
    tags = [
        "no-android",
    ],
    deps = [
        ":selection",
        ":selection_test_util",
        "//seurat/base",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "selection_problem_test",
    size = "small",
    srcs = [
        "selection_problem_test.cc",
    ],
    tags = [
        "no-android",
    ],
    deps = [
        ":selection",
        ":selection_test_util",
        "//seurat/base",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "selection_util_test",
    size = "small",
    srcs = [
        "selection_util_test.cc",
    ],
    tags = [
        "no-android",
    ],
    deps = [
        ":selection",
        "//seurat/base",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "simplifying_solver_test",
    size = "small",
    srcs = [
        "simplifying_solver_test.cc",
    ],
    tags = [
        "no-android",
    ],
    deps = [
        ":selection",
        ":selection_test_util",
        "//seurat/base",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "solver_factory_test",
    size = "small",
    srcs = [
        "solver_factory_test.cc",
    ],
    tags = [
        "no-android",
    ],
    deps = [
        ":selection",
        ":selection_test_util",
        "//seurat/base",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "subgradient_descent_solver_test",
    size = "small",
    srcs = [
        "subgradient_descent_solver_test.cc",
    ],
    tags = [
        "no-android",
    ],
    deps = [
        ":selection",
        ":selection_test_util",
        "//seurat/base",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

test_suite(
    name = "selection_tests",
    tests = [
        ":bisection_solver_test",
        ":cost_calculator_test",
        ":parallel_solver_test",
        ":relaxed_solver_test",
        ":selection_problem_test",
        ":selection_util_test",
        ":simplifying_solver_test",
        ":solver_factory_test",
        ":subgradient_descent_solver_test",
    ],
)
