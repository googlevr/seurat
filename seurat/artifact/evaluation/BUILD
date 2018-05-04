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

# This package includes utilities for evaluating quality of SoftServe Artifacts.

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
    name = "evaluation",
    srcs = [
        "artifact_evaluator.cc",
        "cost_estimator.cc",
        "render_sim.cc",
    ],
    hdrs = [
        "artifact_evaluator.h",
        "cost_estimator.h",
        "render_sim.h",
    ],
    deps = [
        "//seurat/artifact",
        "//seurat/base",
        "//seurat/geometry",
        "//seurat/geometry:raytracer",
        "//seurat/image",
        "//seurat/ingest",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
        "@ion//:iongfx",
        "@ion//:iongfxutils",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

cc_test(
    name = "artifact_evaluator_test",
    size = "small",
    srcs = ["artifact_evaluator_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":evaluation",
        "//seurat/artifact",
        "//seurat/base",
        "//seurat/component",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/ingest",
        "//seurat/ingest:ingest_test_utils",
        "//seurat/mesh:mesh_component_util",
        "//seurat/testing:test_flags",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:iongfx",
        "@ion//:ionmath",
        "@ion//:statetable",
    ],
)

cc_test(
    name = "render_sim_test",
    size = "small",
    srcs = ["render_sim_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":evaluation",
        "//seurat/base",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:iongfx",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "cost_estimator_test",
    size = "small",
    srcs = ["cost_estimator_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":evaluation",
        "//seurat/geometry",
        "//seurat/mesh:mesh_component_util",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

test_suite(
    name = "evaluation_tests",
    tests = [
        ":artifact_evaluator_test",
        ":cost_estimator_test",
        ":render_sim_test",
    ],
)
