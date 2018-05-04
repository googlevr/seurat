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

package(
    default_visibility = ["//visibility:public"],
)

config_setting(
    name = "windows_x86_64",
    values = {"cpu": "x64_windows"},
    visibility = ["//visibility:public"],
)

cc_library(
    name = "embree",
    srcs = glob(
        ["**/*.cpp"],
        exclude = [
            # Exclude cmake subdirectory.
            #
            # It contains a main() function in check_isa.cpp. If it was
            # included in the sources, it would cause all tests to run this
            # main() function and silently exit.
            "common/cmake/**/*",
            # Exclude all tutorials.
            "tutorials/**/*",
            # Exclude all AVX, AVX2 and AVX512 code. We only build the SSE
            # version.
            "**/*avx*.cpp",
            "kernels/bvh/*8*.cpp",
            "kernels/bvh/*16*.cpp",
            # Exclude unused non-portable code.
            "common/sys/library.cpp",
            "common/sys/network.cpp",
            # Exclude unused task systems.
            "common/tasking/taskschedulertbb.cpp",
            "common/tasking/taskschedulerppl.cpp",
        ],
    ),
    hdrs = glob(["include/embree2/*.h"]),
    copts = select({
        "windows_x86_64": ["/EHsc"],
        "//conditions:default": ["-fexceptions"],
    }) + [
        "-w",
        "-DEMBREE_TARGET_SSE2",
        "-D__SSE__",
        "-D__SSE2__",
        # Required to build with Mac.
        "-U__SSE4_1__",
        # Required for Embree to ifdef enable shared code, that would only go
        # into the lowest ISA lib, if multiple ISA libs are built.
        "-DEMBREE_LOWEST_ISA",
    ],
    defines = [
        "TASKING_INTERNAL",
        "EMBREE_STATIC_LIB",
    ],
    # Required because many relative paths in Embree are incorrect, refering
    # to kernels/common instead of common/
    includes = ["common/"],
    linkopts = select({
        ":windows_x86_64": [
            # Required for some Windows API calls in common/sys/alloc.cpp
            "-DEFAULTLIB:advapi32.lib",
        ],
        "//conditions:default": [
        ],
    }),
    textual_hdrs = glob(
        [
            "**/*.h",
            # This .cpp file is included by other .cpp files.
            "kernels/bvh/bvh_intersector1.cpp",
            "kernels/bvh/bvh_intersector_stream.cpp",
            "kernels/bvh/bvh_intersector_hybrid.cpp",
        ],
        exclude = [
            # Exclude all tutorials.
            "tutorials/**/*",
            # Exclude all AVX, AVX2 and AVX512 code. We only build the SSE
            # version.
            "kernels/bvh/*16*.h",
            # Exclude unused non-portable code.
            "common/sys/library.h",
            "common/sys/network.h",
            # Exclude unused task systems.
            "common/tasking/taskschedulertbb.h",
            "common/tasking/taskschedulerppl.h",
        ],
    ),
    deps = [
        ":embree_config",
        ":embree_hash",
    ],
)

#
# When Embree is built with CMake, it generates common/config.h and
# common/hash.h
#
# To make this work with Bazel, versions of those files are inlined here,
# along with genrule/cc_library rules to paste them into files & include them
# as headers.
#

# Generate "common/config.h"
cc_library(
    name = "embree_config",
    srcs = [],
    hdrs = ["common/config.h"],
    visibility = ["//visibility:private"],
)

COMMON_CONFIG_H = """
/* #undef EMBREE_RAY_MASK */
/* #undef EMBREE_STAT_COUNTERS */
/* #undef EMBREE_BACKFACE_CULLING */
#define EMBREE_INTERSECTION_FILTER
#define EMBREE_INTERSECTION_FILTER_RESTORE
/* #undef EMBREE_RETURN_SUBDIV_NORMAL */
#define EMBREE_IGNORE_INVALID_RAYS
#define EMBREE_GEOMETRY_TRIANGLES
#define EMBREE_GEOMETRY_QUADS
/* #undef EMBREE_GEOMETRY_LINES */
/* #undef EMBREE_GEOMETRY_HAIR */
/* #undef EMBREE_GEOMETRY_SUBDIV */
/* #undef EMBREE_GEOMETRY_USER */
#define EMBREE_RAY_PACKETS
/* #undef EMBREE_NATIVE_CURVE_BSPLINE */

#if defined(EMBREE_GEOMETRY_TRIANGLES)
  #define IF_ENABLED_TRIS(x) x
#else
  #define IF_ENABLED_TRIS(x)
#endif

#if defined(EMBREE_GEOMETRY_QUADS)
  #define IF_ENABLED_QUADS(x) x
#else
  #define IF_ENABLED_QUADS(x)
#endif

#if defined(EMBREE_GEOMETRY_LINES)
  #define IF_ENABLED_LINES(x) x
#else
  #define IF_ENABLED_LINES(x)
#endif

#if defined(EMBREE_GEOMETRY_HAIR)
  #define IF_ENABLED_HAIR(x) x
#else
  #define IF_ENABLED_HAIR(x)
#endif

#if defined(EMBREE_GEOMETRY_SUBDIV)
  #define IF_ENABLED_SUBDIV(x) x
#else
  #define IF_ENABLED_SUBDIV(x)
#endif

#if defined(EMBREE_GEOMETRY_USER)
  #define IF_ENABLED_USER(x) x
#else
  #define IF_ENABLED_USER(x)
#endif
"""

genrule(
    name = "generate_config",
    outs = ["common/config.h"],
    cmd = ("echo '%s' > $(OUTS)" % COMMON_CONFIG_H),
    visibility = ["//visibility:private"],
)

# Generate "common/hash.h"
cc_library(
    name = "embree_hash",
    srcs = [],
    hdrs = ["common/hash.h"],
    visibility = ["//visibility:private"],
)

COMMON_HASH_H = """
#define RTCORE_HASH "fbeaa1002fc9c328d14ef06ade213395aff51033"
"""

genrule(
    name = "generate_hash",
    outs = ["common/hash.h"],
    cmd = ("echo '%s' > $(OUTS)" % COMMON_HASH_H),
    visibility = ["//visibility:private"],
)
