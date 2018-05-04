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
    name = "artifact",
    srcs = [
        "atlas_mesh_transform.cc",
        "compress_atlas_transform.cc",
        "compress_tiles_transform.cc",
        "exr_exporter.cc",
        "ice_exporter.cc",
        "mesh_artifact_util.cc",
        "obj_exporter.cc",
        "png_exporter.cc",
        "separate_opaque.cc",
        "sort_atlas_tiles_transform.cc",
    ],
    hdrs = [
        "artifact.h",
        "artifact_processor.h",
        "atlas_mesh_transform.h",
        "compress_atlas_transform.h",
        "compress_tiles_transform.h",
        "exr_exporter.h",
        "ice_exporter.h",
        "mesh_artifact_util.h",
        "obj_exporter.h",
        "png_exporter.h",
        "separate_opaque.h",
        "sort_atlas_tiles_transform.h",
    ],
    deps = [
        "//seurat/base",
        "//seurat/component",
        "//seurat/compressor/rgba",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/image:compression",
        "//seurat/image:exr_io",
        "//seurat/ingest",
        "//seurat/mesh:mesh_component_util",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
        "@ion//:ionbase",
        "@ion//:iongfx",
        "@ion//:ionimage",
        "@ion//:ionmath",
        "@ion//:ionportgfx",
    ],
)

# Artifact test helpers.
cc_library(
    name = "artifact_test_util",
    testonly = 1,
    srcs = [
        "artifact_test_util.cc",
    ],
    hdrs = [
        "artifact_test_util.h",
    ],
    deps = [
        ":artifact",
        "//seurat/base",
        "//seurat/component",
        "//seurat/geometry",
        "//seurat/image",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
        "@ion//:ionbase",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "exporter_test",
    size = "small",
    srcs = ["exporter_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":artifact",
        "//seurat/base",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/mesh:mesh_component_util",
        "//seurat/testing:test_flags",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "artifact_processor_test",
    size = "small",
    srcs = ["artifact_processor_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":artifact",
        "//seurat/base",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/mesh:mesh_component_util",
        "//seurat/testing",
        "//seurat/testing:test_flags",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "compress_atlas_transform_test",
    size = "small",
    srcs = ["compress_atlas_transform_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":artifact",
        ":artifact_test_util",
        "//seurat/base",
        "//seurat/component",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/image:compression",
        "//seurat/image:image_test_utils",
        "//seurat/mesh:mesh_component_util",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:iongfx",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "compress_tiles_transform_test",
    size = "small",
    srcs = ["compress_tiles_transform_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":artifact",
        ":artifact_test_util",
        "//seurat/base",
        "//seurat/component",
        "//seurat/compressor/rgba",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/image:image_test_utils",
        "//seurat/mesh:mesh_component_util",
        "//seurat/testing",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:iongfx",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "atlas_mesh_transform_test",
    size = "small",
    srcs = ["atlas_mesh_transform_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":artifact",
        ":artifact_test_util",
        "//seurat/base",
        "//seurat/component",
        "//seurat/compressor/rgba",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/mesh:mesh_component_util",
        "//seurat/testing",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:iongfx",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "separate_opaque_test",
    size = "small",
    srcs = ["separate_opaque_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":artifact",
        ":artifact_test_util",
        "//seurat/base",
        "//seurat/component",
        "//seurat/compressor/rgba",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/mesh:mesh_component_util",
        "//seurat/testing",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:iongfx",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "mesh_artifact_util_test",
    size = "small",
    srcs = ["mesh_artifact_util_test.cc"],
    tags = [
        "no-android",
    ],
    deps = [
        ":artifact",
        ":artifact_test_util",
        "//seurat/base",
        "//seurat/component",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/mesh:mesh_component_util",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:iongfx",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "sort_atlas_tiles_transform_test",
    size = "small",
    srcs = ["sort_atlas_tiles_transform_test.cc"],
    deps = [
        ":artifact",
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

test_suite(
    name = "artifact_tests",
    tests = [
        ":artifact_processor_test",
        ":atlas_mesh_transform_test",
        ":compress_atlas_transform_test",
        ":compress_tiles_transform_test",
        ":exporter_test",
        ":mesh_artifact_util_test",
        ":separate_opaque_test",
        ":sort_atlas_tiles_transform_test",
    ],
)
