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
    name = "tiler",
    srcs = [
        "build_partition.cc",
        "candidate_tile_generator.cc",
        "geometry_solver_util.cc",
        "partitioner_stages.cc",
        "plane_projection_cost_function.cc",
        "rail_disk_solver.cc",
        "rail_penalty_cost_function.cc",
        "selection_tiler.cc",
        "subdivision.cc",
        "subdivision_util.cc",
        "subset_geometry_solver.cc",
        "tangential_disk_cost_function.cc",
        "tile_resolver.cc",
        "tile_weight_model.cc",
        "tiler.cc",
    ],
    hdrs = [
        "build_partition.h",
        "candidate_tile_generator.h",
        "geometry_model.h",
        "geometry_solver.h",
        "geometry_solver_util.h",
        "partitioner_stage.h",
        "partitioner_stages.h",
        "plane_projection_cost_function.h",
        "point_set.h",
        "rail_disk_solver.h",
        "rail_penalty_cost_function.h",
        "selection_tiler.h",
        "subdivision.h",
        "subdivision_util.h",
        "subset_geometry_solver.h",
        "tangential_disk_cost_function.h",
        "tile.h",
        "tile_resolver.h",
        "tile_weight_model.h",
        "tiler.h",
    ],
    deps = [
        "//seurat/base",
        "//seurat/geometry",
        "//seurat/geometry:kdtree",
        "//seurat/tiler/selection",
        "@ceres_solver//:ceres",
        "@com_google_absl//absl/types:span",
        "@ion//:ionbase",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

cc_library(
    name = "tiler_test_util",
    testonly = 1,
    srcs = [
        "tiler_test_util.cc",
    ],
    hdrs = [
        "tiler_test_util.h",
    ],
    deps = [
        ":tiler",
        "//seurat/base",
        "//seurat/geometry",
        "//seurat/geometry:raytracer",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

cc_test(
    name = "build_partition_test",
    size = "small",
    srcs = ["build_partition_test.cc"],
    deps = [
        ":tiler",
        "//seurat/geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "tiler_test",
    # The full tiler does not work well with very small problems, so
    # tests may take a long time to execute.
    size = "large",
    srcs = ["tiler_test.cc"],
    deps = [
        ":tiler",
        ":tiler_test_util",
        "//seurat/base",
        "//seurat/geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "partitioner_stages_test",
    size = "small",
    srcs = ["partitioner_stages_test.cc"],
    deps = [
        ":tiler",
        ":tiler_test_util",
        "//seurat/base",
        "//seurat/geometry",
        "//third_party/absl/random",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "rail_disk_solver_test",
    size = "small",
    srcs = ["rail_disk_solver_test.cc"],
    deps = [
        ":tiler",
        "//seurat/geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "tile_resolver_test",
    size = "small",
    srcs = ["tile_resolver_test.cc"],
    deps = [
        ":tiler",
        "//seurat/base",
        "//seurat/geometry",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "tile_weight_model_test",
    size = "small",
    srcs = ["tile_weight_model_test.cc"],
    deps = [
        ":tiler",
        "//seurat/base",
        "//seurat/geometry",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "subdivision_test",
    size = "small",
    srcs = ["subdivision_test.cc"],
    deps = [
        ":tiler",
        ":tiler_test_util",
        "//seurat/base",
        "//seurat/geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "subdivision_util_test",
    size = "small",
    srcs = ["subdivision_util_test.cc"],
    deps = [
        ":tiler",
        "//seurat/base",
        "//seurat/geometry",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "subset_geometry_solver_test",
    size = "small",
    srcs = ["subset_geometry_solver_test.cc"],
    deps = [
        ":tiler",
        "//seurat/base",
        "//seurat/geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

test_suite(
    name = "tiler_tests",
    tests = [
        ":build_partition_test",
        ":partitioner_stages_test",
        ":rail_disk_solver_test",
        ":subdivision_test",
        ":subdivision_util_test",
        ":subset_geometry_solver_test",
        ":tile_resolver_test",
        ":tile_weight_model_test",
        ":tiler_test",
    ],
)
