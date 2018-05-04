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
    name = "framework",
    srcs = [
        "frame.cc",
        "frame_generator.cc",
        "frame_sorter.cc",
        "radiance_accumulator.cc",
        "rasterizer.cc",
        "ray_bundle.cc",
        "ray_classifier.cc",
        "silhouette.cc",
        "silhouette_accumulator.cc",
        "texture_sizer.cc",
    ],
    hdrs = [
        "frame.h",
        "frame_generator.h",
        "frame_sorter.h",
        "radiance_accumulator.h",
        "rasterizer.h",
        "ray_bundle.h",
        "ray_classifier.h",
        "sample_accumulator.h",
        "sample_accumulator_set.h",
        "silhouette.h",
        "silhouette_accumulator.h",
        "texture_sizer.h",
    ],
    deps = [
        "//seurat/base",
        "//seurat/geometry",
        "//seurat/geometry:kdtree",
        "//seurat/geometry:raytracer",
        "//seurat/image",
        "//seurat/ingest",
        "//seurat/mesh",
        "//seurat/mesh:mesh_component_util",
        "//seurat/tiler",
        "@com_google_absl//absl/types:span",
        "@ion//:iongfx",
        "@ion//:iongfxutils",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

cc_test(
    name = "frame_test",
    size = "small",
    srcs = ["frame_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":framework",
        "//seurat/base",
        "//seurat/image",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "frame_generator_test",
    # This test processes a tiny scene into frames, including running the
    # Tiler.
    size = "enormous",
    srcs = ["frame_generator_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":framework",
        "//seurat/base",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/ingest",
        "//seurat/ingest:ingest_test_utils",
        "//seurat/mesh",
        "//seurat/mesh:mesh_component_util",
        "//seurat/testing",
        "//seurat/testing:test_flags",
        "//seurat/tiler",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "ray_classifier_test",
    size = "small",
    srcs = ["ray_classifier_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":framework",
        "//seurat/base",
        "//seurat/image",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "frame_sorter_test",
    size = "small",
    srcs = ["frame_sorter_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":framework",
        "//seurat/base",
        "//seurat/geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "silhouette_test",
    size = "small",
    srcs = ["silhouette_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":framework",
        "//seurat/base",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "texture_sizer_test",
    size = "small",
    srcs = ["texture_sizer_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":framework",
        "//seurat/base",
        "//seurat/image",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

test_suite(
    name = "framework_tests",
    tests = [
        ":frame_generator_test",
        ":frame_sorter_test",
        ":frame_test",
        ":ray_classifier_test",
        ":silhouette_test",
        ":texture_sizer_test",
    ],
)
