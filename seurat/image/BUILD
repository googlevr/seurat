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

config_setting(
    name = "linux_x86_64",
    values = {"cpu": "k8"},
)

config_setting(
    name = "windows_x86_64",
    values = {"cpu": "x64_windows"},
)

cc_library(
    name = "image",
    srcs = [
        "atlaser_util.cc",
        "color_processor.cc",
        "fixed_width_atlaser.cc",
        "image_util.cc",
        "inpainting.cc",
        "ion_image_reader.cc",
        "ion_image_writer.cc",
        "ldi_util.cc",
        "nearly_square_atlaser.cc",
        "texture_atlas.cc",
    ],
    hdrs = [
        "atlaser.h",
        "atlaser_util.h",
        "color_processor.h",
        "filter.h",
        "fixed_width_atlaser.h",
        "image.h",
        "image_reader.h",
        "image_util.h",
        "inpainting.h",
        "ion_image_reader.h",
        "ion_image_writer.h",
        "ldi.h",
        "ldi_util.h",
        "low_discrepancy_sampling.h",
        "nearly_square_atlaser.h",
        "texture_atlas.h",
    ],
    deps = [
        "//seurat/base",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
        "@ion//:ionbase",
        "@ion//:iongfx",
        "@ion//:ionimage",
        "@ion//:ionmath",
        "@ion//:ionportgfx",
    ],
)

cc_library(
    name = "compression",
    srcs = [
        "rgbaui8_codec.cc",
    ],
    hdrs = [
        "codec.h",
        "rgbaui8_codec.h",
    ],
    deps = [
        ":image",
        "//seurat/base",
        "@ion//:ionbase",
        "@ion//:iongfx",
        "@ion//:ionimage",
        "@ion//:ionmath",
    ],
)

EXCEPTIONS_COPTS = select({
    "linux_x86_64": ["-fexceptions"],
    "windows_x86_64": ["/EHsc"],
})

# The EXR I/O code is isolated in its own target since it needs to be built with
# exceptions.  This target is depended on by :image, so it cannot have circular
# dependencies back to :image.  In fact, it is written to have as few
# dependencies as possible.
cc_library(
    name = "exr_io",
    srcs = [
        "deep_exr_image_reader.cc",
        "deep_exr_image_writer.cc",
        "exr_image_reader.cc",
        "exr_image_writer.cc",
        "exr_util.cc",
        "image_reader.h",
        "image_writer.h",
    ],
    hdrs = [
        "deep_exr_image_reader.h",
        "deep_exr_image_writer.h",
        "exr_image_reader.h",
        "exr_image_writer.h",
        "exr_util.h",
    ],
    copts = EXCEPTIONS_COPTS,
    deps = [
        "//seurat/base",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
        "@openexr//:half",
        "@openexr//:iex",
        "@openexr//:ilm_imf",
        "@openexr//:imath",
    ],
)

cc_library(
    name = "image_test_utils",
    testonly = 1,
    srcs = [
        "image_test_utils.cc",
    ],
    hdrs = [
        "image_test_utils.h",
        "ldi_test_utils.h",
    ],
    deps = [
        ":image",
        "//seurat/api:image_portable_proto",
        "//seurat/base",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "filter_test",
    size = "small",
    srcs = ["filter_test.cc"],
    deps = [
        ":image",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "image_test",
    size = "small",
    srcs = ["image_test.cc"],
    deps = [
        ":image",
        "//seurat/base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "ion_io_test",
    size = "small",
    srcs = ["ion_io_test.cc"],
    data = [
        ":testdata/color.png",
    ],
    tags = ["no-android"],
    deps = [
        ":image",
        "//seurat/base",
        "//seurat/testing",
        "//seurat/testing:test_flags",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:iongfx",
        "@ion//:ionimage",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "image_util_test",
    size = "small",
    srcs = ["image_util_test.cc"],
    deps = [
        ":image",
        "//seurat/base",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "ldi_test",
    size = "small",
    srcs = ["ldi_test.cc"],
    deps = [
        ":image",
        ":image_test_utils",
        "//seurat/base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "ldi_util_test",
    size = "small",
    srcs = ["ldi_util_test.cc"],
    deps = [
        ":image",
        "//seurat/base",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "low_discrepancy_sampling_test",
    size = "small",
    srcs = ["low_discrepancy_sampling_test.cc"],
    deps = [
        ":image",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "inpainting_test",
    size = "small",
    srcs = ["inpainting_test.cc"],
    deps = [
        ":image",
        "//seurat/base",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "nearly_square_atlaser_test",
    size = "small",
    srcs = [
        "nearly_square_atlaser_test.cc",
    ],
    deps = [
        ":image",
        "//seurat/base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "fixed_width_atlaser_test",
    size = "small",
    srcs = [
        "fixed_width_atlaser_test.cc",
    ],
    deps = [
        ":image",
        "//seurat/base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "texture_atlas_test",
    size = "small",
    srcs = ["texture_atlas_test.cc"],
    deps = [
        ":image",
        "//seurat/base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "deep_exr_io_test",
    size = "small",
    srcs = [
        "deep_exr_io_test.cc",
    ],
    deps = [
        ":exr_io",
        ":image",
        ":image_test_utils",
        "//seurat/base",
        "//seurat/testing",
        "//seurat/testing:test_flags",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "exr_io_test",
    size = "small",
    srcs = [
        "exr_io_test.cc",
    ],
    data = [
        ":testdata/color.exr",
        ":testdata/depth.exr",
    ],
    deps = [
        ":exr_io",
        ":image",
        ":image_test_utils",
        "//seurat/base",
        "//seurat/testing",
        "//seurat/testing:test_flags",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "color_processor_test",
    size = "small",
    srcs = ["color_processor_test.cc"],
    deps = [
        ":image",
        "//seurat/base",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "rgbaui8_codec_test",
    size = "small",
    srcs = ["rgbaui8_codec_test.cc"],
    deps = [
        ":compression",
        ":image",
        "//seurat/base",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

cc_test(
    name = "atlaser_util_test",
    size = "small",
    srcs = ["atlaser_util_test.cc"],
    deps = [
        ":image",
        "//seurat/base",
        "//seurat/testing",
        "@com_google_googletest//:gtest",
        "@com_google_googletest//:gtest_main",
        "@ion//:ionmath",
    ],
)

test_suite(
    name = "image_tests",
    tests = [
        ":atlaser_util_test",
        ":color_processor_test",
        ":deep_exr_io_test",
        ":exr_io_test",
        ":filter_test",
        ":fixed_width_atlaser_test",
        ":image_test",
        ":image_util_test",
        ":inpainting_test",
        ":ion_io_test",
        ":ldi_test",
        ":ldi_util_test",
        ":low_discrepancy_sampling_test",
        ":nearly_square_atlaser_test",
        ":rgbaui8_codec_test",
        ":texture_atlas_test",
    ],
)
