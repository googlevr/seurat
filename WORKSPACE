workspace(name = "com_google_seurat")

local_repository(
    name = "com_google_absl",
    path = "third_party/ion/third_party/absl",
)

# CCTZ (Time-zone framework).
http_archive(
    name = "com_googlesource_code_cctz",
    strip_prefix = "cctz-2.1",
    urls = ["https://github.com/google/cctz/archive/v2.1.zip"],
)

http_archive(
    name = "com_google_protobuf",
    strip_prefix = "protobuf-3.4.1",
    urls = ["https://github.com/google/protobuf/archive/v3.4.1.zip"],
)

http_archive(
    name = "com_google_protobuf_cc",
    strip_prefix = "protobuf-3.4.1",
    urls = ["https://github.com/google/protobuf/archive/v3.4.1.zip"],
)

# GoogleTest/GoogleMock framework. Used by most unit-tests.
http_archive(
    name = "com_google_googletest",
    strip_prefix = "googletest-master",
    urls = ["https://github.com/google/googletest/archive/master.zip"],
)

# Khronos GL headers.
new_http_archive(
    name = "khronos_opengl_registry",
    strip_prefix = "OpenGL-Registry-master",
    build_file = "third_party/gl.BUILD",
    urls = ["https://github.com/KhronosGroup/OpenGL-Registry/archive/master.zip"],
)

new_local_repository(
    name = "ion",
    build_file = "third_party/ion/BUILD",
    path = "third_party/ion/",
)

new_http_archive(
    name = "eigen",
    build_file = "third_party/eigen.BUILD",
    strip_prefix = "eigen-eigen-5a0156e40feb",
    urls = ["https://bitbucket.org/eigen/eigen/get/3.3.4.zip"],
)

new_http_archive(
    name = "ceres_solver",
    build_file = "third_party/ceres.BUILD",
    strip_prefix = "ceres-solver-1.13.0",
    urls = ["https://github.com/ceres-solver/ceres-solver/archive/1.13.0.zip"],
)

new_http_archive(
    name = "nanoflann",
    build_file = "third_party/nanoflann.BUILD",
    strip_prefix = "nanoflann-1.2.3",
    urls = ["https://github.com/jlblancoc/nanoflann/archive/v1.2.3.zip"],
)

new_http_archive(
    name = "tbb",
    build_file = "third_party/tbb.BUILD",
    strip_prefix = "tbb-2017_U1",
    urls = ["https://github.com/01org/tbb/archive/2017_U1.zip"],
)

new_http_archive(
    name = "embree",
    build_file = "third_party/embree.BUILD",
    strip_prefix = "embree-2.16.5",
    urls = ["https://github.com/embree/embree/archive/v2.16.5.zip"],
)

http_archive(
    name = "gflags",
    strip_prefix = "gflags-6d1c363dde8a4277f7f363446a82b2188fd94608",
    url = "https://github.com/gflags/gflags/archive/6d1c363dde8a4277f7f363446a82b2188fd94608.zip",
)

new_http_archive(
    name = "openexr",
    build_file = "third_party/openexr.BUILD",
    strip_prefix = "openexr-2.2.0",
    urls = ["https://github.com/openexr/openexr/archive/v2.2.0.zip"],
)

new_http_archive(
    name = "glfw",
    build_file = "third_party/glfw.BUILD",
    strip_prefix = "glfw-3.2.1",
    urls = ["https://github.com/glfw/glfw/archive/3.2.1.zip"],
)
