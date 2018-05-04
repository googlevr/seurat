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

CERES_DEFAULT_DEFINES = [
    # Disabling c++11 support improves portability between different
    # compilers (due to inconstent "support" for various C++11 features).
    #
    # "CERES_USE_CXX11",

    "CERES_NO_SUITESPARSE",

    # Don't enable support for CXSparse.
    "CERES_NO_CXSPARSE",

    # Unordered map is in the std namespace.
    "CERES_STD_UNORDERED_MAP",

    # Do not use GPL Code in Eigen, including Eigen's Sparse Cholesky
    # module.
    "EIGEN_MPL2_ONLY",
    "CERES_NO_EIGENSPARSE",

    # Do not call lapack directly (this ensures, FORTRAN LAPACK is not
    # needed).
    "CERES_NO_LAPACK",
]

CERES_THREADING_DEFINES = [
    "CERES_NO_THREADS",
    "CERES_HAVE_PTHREAD",
    "CERES_HAVE_RWLOCK",
]

CERES_DEFINES = CERES_DEFAULT_DEFINES + CERES_THREADING_DEFINES

cc_library(
    name = "ceres_config",
    srcs = [],
    hdrs = glob([
        "config/ceres/internal/*.h",
    ]),
    defines = CERES_DEFINES,
    includes = ["config"],
)

cc_library(
    name = "ceres_glog",
    srcs = [
        "include/ceres/internal/disable_warnings.h",
        "include/ceres/internal/port.h",
        "include/ceres/internal/reenable_warnings.h",
        "internal/ceres/miniglog/glog/logging.cc",
    ],
    hdrs = ["internal/ceres/miniglog/glog/logging.h"],
    includes = ["include"],
    strip_include_prefix = "internal/ceres/miniglog",
    deps = [
        ":ceres_config",
    ],
)

cc_library(
    name = "ceres",
    srcs = glob(
        [
            "internal/ceres/*.cc",
            "internal/ceres/*.h",
            "internal/ceres/generated/*.cc",
            "internal/ceres/generated/*.h",
        ],
        exclude = [
            "internal/ceres/*_test.cc",
            "internal/ceres/*_test_utils.cc",
            "internal/ceres/*_test_utils.h",
            "internal/ceres/gmock*",
            "internal/ceres/gtest*",
            "internal/ceres/test_util.*",
        ],
    ),
    hdrs = glob([
        "include/ceres/*.h",
        "include/ceres/internal/*.h",
    ]),
    copts = [
        # Disable verbose logging.
        "-DMAX_LOG_LEVEL=0",
    ],
    defines = CERES_DEFINES,
    includes = [
        "include",
        "internal",
    ],
    linkstatic = 1,
    deps = [
        ":ceres_config",
        ":ceres_glog",
        "@eigen//:eigen",
    ],
)
