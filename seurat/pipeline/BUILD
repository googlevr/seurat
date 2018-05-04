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
    "portable_proto_library",
)

package_group(
    name = "friends",
    packages = [
        "//seurat/...",
    ],
)

package(default_visibility = [":friends"])

distribs(["client"])

proto_library(
    name = "flags_proto",
    srcs = ["flags.proto"],
)

cc_proto_library(
    name = "flags_portable_proto",
    deps = ["flags_proto"],
)

pipeline_extra_deps = []

cc_library(
    name = "pipeline",
    srcs = [
        "pipeline.cc",
    ],
    hdrs = [
        "pipeline.h",
    ],
    deps = [
        ":flags_portable_proto",
        "@ion//:ionbase",
        "@ion//:ionmath",
        "@com_google_absl//absl/strings",
        "//seurat/artifact",
        "//seurat/artifact/evaluation",
        "//seurat/baker/diffuse",
        "//seurat/baker/framework",
        "//seurat/base",
        "//seurat/component",
        "//seurat/compressor/rgba",
        "//seurat/geometry",
        "//seurat/image",
        "//seurat/image:compression",
        "//seurat/ingest",
        "//seurat/mesh",
        "//seurat/mesh:mesh_component_util",
        "//seurat/tiler",
    ] + pipeline_extra_deps,
)

seurat_extra_deps = [
    "@gflags//:gflags",
    "@ion//:ion_jsoncpp",
]

cc_binary(
    name = "seurat",
    srcs = ["seurat.cc"],
    deps = [
        ":pipeline",
        "@ion//:ionbase",
        "@ion//:ionport",
        "@com_google_absl//absl/strings",
        "@openexr//:ilm_imf",
        "//seurat/artifact",
        "//seurat/base",
    ] + seurat_extra_deps,
)
