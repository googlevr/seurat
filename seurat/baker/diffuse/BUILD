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
    name = "diffuse",
    srcs = [
        "diffuse_baker.cc",
        "diffuse_pipeline.cc",
    ],
    hdrs = [
        "diffuse_baker.h",
        "diffuse_pipeline.h",
    ],
    deps = [
        "//seurat/artifact",
        "//seurat/baker/framework",
        "//seurat/base",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/mesh",
        "//seurat/mesh:mesh_component_util",
        "@com_google_absl//absl/types:span",
        "@ion//:iongfx",
        "@ion//:iongfxutils",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

# This executes an end-to-end pipeline over a tiny scene.
cc_test(
    name = "diffuse_pipeline_test",
    size = "large",
    srcs = ["diffuse_pipeline_test.cc"],
    tags = [
        "no-android",
        "no-msvc",
    ],
    deps = [
        ":diffuse",
        "//seurat/artifact",
        "//seurat/baker/framework",
        "//seurat/base",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/ingest",
        "//seurat/ingest:ingest_test_utils",
        "//seurat/mesh",
        "//seurat/mesh:mesh_component_util",
        "//seurat/tiler",
        "@com_google_absl//absl/memory",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

test_suite(
    name = "diffuse_tests",
    tests = [
        ":diffuse_pipeline_test",
    ],
)
