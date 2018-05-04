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

# Description:
#   Bazel build file for the ILM OpenEXR libraries.

package(
    default_visibility = ["//visibility:public"],
    features = [
        "-layering_check",
        "-parse_headers",
    ],
)

licenses(["notice"])  # BSD

exports_files(["LICENSE"])

config_setting(
    name = "linux_x86_64",
    values = {"cpu": "k8"},
    visibility = ["//visibility:public"],
)

config_setting(
    name = "windows_x86_64",
    values = {"cpu": "x64_windows"},
    visibility = ["//visibility:public"],
)

DEFINES = [
    "DISABLE_GOOGLE_GLOBAL_USING_DECLARATIONS",
] + select({
    ":linux_x86_64": [
        "OS_LINUX=OS_LINUX",
        "ARCH_K8",
    ],
    # The following match the build environment we've used to build Seurat
    # libraries for Windows.
    ":windows_x86_64": [
        "NOGDI",  # Don't pollute with GDI macros in windows.h.
        "NOMINMAX",  # Don't define min/max macros in windows.h.
        "OS_WINDOWS=OS_WINDOWS",
        "PRAGMA_SUPPORTED",
        "WIN32",
        "WINVER=0x0600",
        "_CRT_SECURE_NO_DEPRECATE",
        "_WIN32",
        "_WIN32_WINNT=0x0600",
        "_WINDOWS",

        # Defaults for other headers (along with OS_WINDOWS).
        "COMPILER_MSVC",

        # Use math constants (M_PI, etc.) from the math library
        "_USE_MATH_DEFINES",

        # Allows 'Foo&&' (e.g., move constructors).
        "COMPILER_HAS_RVALUEREF",

        # Doublespeak for "don't bloat namespace with incompatible winsock
        # versions that I didn't include".
        # http://msdn.microsoft.com/en-us/library/windows/desktop/ms737629.aspx
        "WIN32_LEAN_AND_MEAN=1",

        "PLATFORM_WINDOWS",
    ],
})

# Options for everything below
all_options = [
    "-D_THREAD_SAFE",
] + select({
    ":windows_x86_64": ["/EHsc"],
    ":linux_x86_64": [
        "-fexceptions",
        "-fno-common",
        "-Wno-unused-variable",
        "-Wno-unused-result",
        "-Wno-string-conversion",
        "-Wno-switch",
        "-Wno-tautological-compare",
	"-Wno-deprecated-register",
        # IlmImf/ImfHuf.cpp:111 wants to put a 21520 byte object on the stack.
        "$(STACK_FRAME_UNLIMITED)",
    ],
})

# Generates the config headers. Note we don't use the autoconf, but the
# libraries still require a header.
genrule(
    name = "configure_ilmbase",
    outs = ["IlmBase/config/IlmBaseConfig.h"],
    cmd = "touch $@",
)

# Wraps the config header auto generation, then injects the desired defines and
# include folder configuration.
cc_library(
    name = "ilmbase_config",
    hdrs = [":configure_ilmbase"],
    defines = DEFINES,
    includes = ["IlmBase/config/"],
)

################################################################################
# Build the half library and its support components
# ILM's libhalf is considered the reference implementation.
half_path = "IlmBase/Half/"

cc_binary(
    name = "toFloat",
    srcs = [half_path + "toFloat.cpp"],
    copts = all_options,
    defines = DEFINES,
)

cc_binary(
    name = "eLut",
    srcs = [half_path + "eLut.cpp"],
    copts = all_options,
    defines = DEFINES,
)

genrule(
    name = "toFloatAutoGen",
    outs = [half_path + "toFloat.h"],
    # Runs the tool from its built location, and captures to the single output
    # of this rule.
    cmd = "$(location :toFloat) > $@",
    tools = [":toFloat"],
)

genrule(
    name = "eLutAutoGen",
    outs = [half_path + "eLut.h"],
    # Runs the tool from its built location, and captures to the single output
    # of this rule.
    cmd = "$(location :eLut) > $@",
    tools = [":eLut"],
)

cc_library(
    name = "half",
    srcs = glob([half_path + "half.cpp"]),
    copts = all_options,
    defines = DEFINES,
    includes = [half_path],
    textual_hdrs = glob([
        half_path + "half.h",
        half_path + "halfExport.h",
        half_path + "halfFunction.h",
        half_path + "halfLimits.h",
    ]) + [
        ":toFloatAutoGen",
        ":eLutAutoGen",
    ],
    deps = [":ilmbase_config"],
)

################################################################################
# Build the Iex library
iex_path = "IlmBase/Iex/"

cc_library(
    name = "iex",
    srcs = glob([iex_path + "*.cpp"]),
    hdrs = glob([iex_path + "*.h"]),
    copts = all_options,
    defines = DEFINES,
    includes = [iex_path],
    deps = [":ilmbase_config"],
)

################################################################################
# Build the IexMath library
iex_math_path = "IlmBase/IexMath/"

cc_library(
    name = "iex_math",
    srcs = glob([iex_math_path + "*.cpp"]),
    hdrs = glob([iex_math_path + "*.h"]),
    copts = all_options,
    includes = [iex_math_path],
    deps = [
        ":half",
        ":iex",
    ],
)

################################################################################
# Build the IlmThread library and its tests components
ilm_thread_path = "IlmBase/IlmThread/"

ilm_thread_win_srcs_pattern = ilm_thread_path + "*Win32.cpp"

ilm_thread_win_hdrs_pattern = ilm_thread_path + "*Win32.h"

cc_library(
    name = "ilm_thread",
    srcs = glob(
        [ilm_thread_path + "*.cpp"],
        exclude = [ilm_thread_win_srcs_pattern],
    ) + select({
        ":windows_x86_64": glob([ilm_thread_win_srcs_pattern]),
        ":linux_x86_64": [],
    }),
    hdrs = glob(
        [ilm_thread_path + "*.h"],
        exclude = [ilm_thread_win_hdrs_pattern],
    ) + select({
        ":windows_x86_64": glob([ilm_thread_win_hdrs_pattern]),
        ":linux_x86_64": [],
    }),
    copts = all_options,
    defines = DEFINES,
    includes = [
        ilm_thread_path,
    ],
    deps = [
        ":iex",
        ":ilmbase_config",
    ],
)

################################################################################
# Build the Imath library
imath_path = "IlmBase/Imath/"

cc_library(
    name = "imath",
    srcs = glob([imath_path + "*.cpp"]),
    hdrs = glob([imath_path + "*.h"]),
    copts = all_options,
    defines = DEFINES,
    includes = [
        imath_path,
    ],
    deps = [
        ":half",
        ":iex",
        ":iex_math",
        ":ilmbase_config",
    ],
)

################################################################################
# Build the IlmImf library
ilm_imf_path = "OpenEXR/IlmImf/"

ilm_imf_util_path = "OpenEXR/IlmImfUtil/"

# Generates the config headers. Note we don't use the autoconf, but the
# libraries still require a header.
genrule(
    name = "configure_openexr",
    outs = ["OpenEXR/config/OpenEXRConfig.h"],
    cmd = "touch $@",
)

# Wraps the config auto generation with the include folder configuration.
cc_library(
    name = "openexr_config",
    hdrs = [":configure_openexr"],
    defines = DEFINES,
    includes = ["OpenEXR/config/"],
)

cc_binary(
    name = "dwaLookups",
    srcs = [
        ilm_imf_path + "dwaLookups.cpp",
    ] +
    # Note we include various headers via glob since we don't have cc_library
    # targets for some of the sub-libraries of OpenEXR.
    glob([
        ilm_imf_util_path + "**/*.h",
        ilm_imf_path + "*.h",
    ]),
    copts = all_options,
    defines = DEFINES,
    includes = [
        imath_path,
        iex_math_path,
        ilm_imf_path,
        ilm_imf_util_path,
    ],
    deps = [
        ":half",
        ":iex",
        ":ilm_thread",
        ":imath",
        ":openexr_config",
    ],
)

cc_binary(
    name = "b44ExpLogTable",
    srcs = [
        ilm_imf_path + "b44ExpLogTable.cpp",
    ] + glob(["**/*.h"]),
    copts = all_options,
    defines = DEFINES,
    includes = [
        half_path,
        imath_path,
        iex_math_path,
        iex_path,
        ilm_imf_path,
        ilm_imf_util_path,
        ilm_thread_path,
    ],
    deps = [
        ":half",
        ":iex",
        ":ilm_thread",
        ":openexr_config",
    ],
)

genrule(
    name = "dwaLookupsAutoGen",
    outs = [ilm_imf_path + "dwaLookups.h"],
    # Runs the tool from its built location, and captures to the single output
    # of this rule.
    cmd = "$(location :dwaLookups) > $@",
    tools = [":dwaLookups"],
)

genrule(
    name = "b44ExpLogTableAutoGen",
    outs = [ilm_imf_path + "b44ExpLogTable.h"],
    # Runs the tool from its built location, and captures to the single output
    # of this rule.
    cmd = "$(location :b44ExpLogTable) > $@",
    tools = [":b44ExpLogTable"],
)

cc_library(
    name = "ilm_imf",
    srcs = glob([ilm_imf_path + "Imf*.cpp"]) + [
        ":dwaLookupsAutoGen",
        ":b44ExpLogTableAutoGen",
    ],
    hdrs = glob([ilm_imf_path + "Imf*.h"]),
    copts = all_options,
    defines = DEFINES,
    includes = [ilm_imf_path],
    deps = [
        ":half",
        ":iex",
        ":iex_math",
        ":ilm_thread",
        ":imath",
        ":openexr_config",
        "@ion//:ion_external_libs",
    ],
)
