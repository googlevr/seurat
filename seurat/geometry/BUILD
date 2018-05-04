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

config_setting(
    name = "windows_x86_64",
    values = {"cpu": "x64_windows"},
    visibility = ["//visibility:public"],
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
    name = "geometry",
    srcs = [
        "binning_point_cloud_builder.cc",
        "convex_hull2d.cc",
        "cube_face.cc",
        "fibonacci_sphere.cc",
        "mesh.cc",
        "mesh_obj_io.cc",
        "mesh_util.cc",
        "polygon.cc",
        "quad_mesh.cc",
        "quad_mesh_util.cc",
        "ray_sphere_intersection.cc",
    ],
    hdrs = [
        "bilinear_interpolator.h",
        "binning_point_cloud_builder.h",
        "convex_hull2d.h",
        "cube_face.h",
        "fibonacci_sphere.h",
        "mesh.h",
        "mesh_obj_io.h",
        "mesh_util.h",
        "plane.h",
        "point_cloud_builder.h",
        "polygon.h",
        "quad.h",
        "quad_mesh.h",
        "quad_mesh_util.h",
        "ray_sphere_intersection.h",
        "triangle.h",
    ],
    deps = [
        "//seurat/base",
        "//seurat/image",
        "@com_google_absl//absl/types:span",
        "@ion//:ionbase",
        "@ion//:ionmath",
    ],
)

# This does not build for android, so use a separate build-target to avoid
# pulling it into the viewer.
cc_library(
    name = "raytracer",
    srcs = [
        "raytracer.cc",
    ],
    hdrs = [
        "raytracer.h",
    ],
    deps = [
        ":geometry",
        "//seurat/base",
        "//seurat/image",
        "@embree",
        "@ion//:ionmath",
    ],
)

cc_library(
    name = "kdtree",
    srcs = [
        "kdtree.cc",
    ],
    hdrs = [
        "kdtree.h",
    ],
    copts = select({
        "windows_x86_64": ["/EHsc"],
        # Required for nanoflann on Linux.
        "//conditions:default": ["-fexceptions"],
    }),
    deps = [
        ":geometry",
        "//seurat/base",
        "@com_google_absl//absl/types:span",
        "@ion//:ionbase",
        "@ion//:ionmath",
        "@nanoflann",
    ],
)

cc_test(
    name = "convex_hull2d_test",
    size = "small",
    srcs = ["convex_hull2d_test.cc"],
    deps = [
        ":geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "polygon_test",
    size = "small",
    srcs = ["polygon_test.cc"],
    deps = [
        ":geometry",
        "//seurat/base",
        "//seurat/image",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionbase",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "cube_face_test",
    size = "small",
    srcs = ["cube_face_test.cc"],
    deps = [
        ":geometry",
        "//seurat/base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionbase",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "mesh_test",
    size = "small",
    srcs = ["mesh_test.cc"],
    deps = [
        ":geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionbase",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "plane_test",
    size = "small",
    srcs = ["plane_test.cc"],
    deps = [
        ":geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionbase",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "fibonacci_sphere_test",
    size = "small",
    srcs = ["fibonacci_sphere_test.cc"],
    deps = [
        ":geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionbase",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "mesh_util_test",
    size = "small",
    srcs = ["mesh_util_test.cc"],
    deps = [
        ":geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionbase",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "mesh_obj_io_test",
    size = "small",
    srcs = ["mesh_obj_io_test.cc"],
    deps = [
        ":geometry",
        "//seurat/base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionbase",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "kdtree_test",
    size = "small",
    srcs = ["kdtree_test.cc"],
    deps = [
        ":geometry",
        ":kdtree",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "triangle_test",
    size = "small",
    srcs = ["triangle_test.cc"],
    deps = [
        ":geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionbase",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "ray_sphere_intersection_test",
    size = "small",
    srcs = ["ray_sphere_intersection_test.cc"],
    deps = [
        ":geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "raytracer_test",
    size = "small",
    srcs = ["raytracer_test.cc"],
    tags = [
        # Test is flaky on msvc.
        "no-msvc",
    ],
    deps = [
        ":geometry",
        ":raytracer",
        "//seurat/base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "binning_point_cloud_builder_test",
    size = "medium",
    srcs = ["binning_point_cloud_builder_test.cc"],
    deps = [
        ":geometry",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

cc_test(
    name = "bilinear_interpolator_test",
    size = "medium",
    srcs = ["bilinear_interpolator_test.cc"],
    deps = [
        ":geometry",
        "//seurat/base",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
        "@ion//:ionport",
    ],
)

cc_test(
    name = "quad_mesh_test",
    size = "small",
    srcs = ["quad_mesh_test.cc"],
    deps = [
        ":geometry",
        "//seurat/base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

test_suite(
    name = "geometry_tests",
    tests = [
        ":bilinear_interpolator_test",
        ":binning_point_cloud_builder_test",
        ":convex_hull2d_test",
        ":cube_face_test",
        ":fibonacci_sphere_test",
        ":kdtree_test",
        ":mesh_obj_io_test",
        ":mesh_test",
        ":mesh_util_test",
        ":plane_test",
        ":polygon_test",
        ":quad_mesh_test",
        ":ray_sphere_intersection_test",
        ":raytracer_test",
        ":triangle_test",
    ],
)
