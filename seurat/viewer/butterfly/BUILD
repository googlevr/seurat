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

package_group(
    name = "friends",
    packages = [
        "//seurat/...",
    ],
)

package(
    default_visibility = [":friends"],
)

config_setting(
    name = "linux",
    constraint_values = [
        "@bazel_tools//platforms:linux",
    ],
)

LINKOPTS = select({
    ":linux": [
        "-ldl",
        "-lX11",
        "-lXrandr",
        "-lXi",
        "-lGL",
        "-lXcursor",
        "-lXinerama",
    ],
    "//conditions:default": [],
})

EXTRA_LIBRARY_DEPS = [
    "@khronos_opengl_registry//:gl_headers",
]

cc_binary(
    name = "butterfly",
    srcs = [
        "butterfly.cc",
        "viewer_camera.cc",
        "viewer_camera.h",
    ],
    copts = [
        "-DGL_GLEXT_PROTOTYPES",
        "-DGLFW_INCLUDE_GLEXT",
    ],
    linkopts = LINKOPTS,
    deps = [
        "@ion//:graphicsmanager",
        "@ion//:iongfx",
        "@ion//:tracinghelper",
        "@ion//:iongfxutils",
        "@ion//:ionmath",
        "@ion//:ionport",
        "@ion//:ionportgfx",
        "@ion//:ionremote",
        "@glfw//:glfw",
        "@com_google_absl//absl/strings",
        "//seurat/base",
        "//seurat/component",
        "//seurat/viewer",
    ] + EXTRA_LIBRARY_DEPS,
)
