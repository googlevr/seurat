# How to Build

## Linux

1. Download Ion into third_party/ion and build via their instructions:
```
~/Seurat$ git clone https://github.com/google/Ion third_party/ion
~/Seurat$ cd third_party/ion
~/Seurat/third_party/ion$ git submodule init
~/Seurat/third_party/ion$ git submodule update
~/Seurat/third_party/ion$ cd ion
~/Seurat/third_party/ion$ ./build.sh -c opt
~/Seurat/third_party/ion$ ./build.sh -c opt external/harfbuzz.gyp
```

2. Build the pipeline executable with:
```
~/Seurat$ bazel --bazelrc tools/bazel_linux.rc build -c opt seurat/pipeline:seurat
```
   Note that "opt" builds can be 100x faster and are necessary for the
   pipeline to complete in reasonable time.

3. Build the viewer executable with
```
~/Seurat$ bazel --bazelrc tools/bazel_linux.rc build -c opt seurat/viewer/butterfly:butterfly
```

### Supported Toolchain

#### Bazel 0.9.0

The build is verified to work with Bazel release 0.9.0. Other versions,
including newer ones, might not work.

#### Clang 3.9

If building on Ubuntu 14.04 LTS, the default of clang-3.4 will not work.

clang-3.9 is known to work.

To use clang-3.9, run the following:

```
$ sudo apt-get install clang-3.9
$ sudo update-alternatives --install /usr/bin/clang clang /usr/bin/clang-3.9 100
$ sudo update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-3.9 100

```

Verify that the correct version is installed by running
```
$ clang --version
```
and
```
$ clang++ --version
```

## Windows

1. Download Ion into third\_party\\ion and build via ion\\ion\\build.bat -c opt_x64
```
c:\Seurat>git clone https://github.com/google/Ion third\_party\\ion
c:\Seurat>cd third\_party\\ion
c:\Seurat\third\_party\\ion>git submodule init
c:\Seurat\third\_party\\ion>git submodule update
c:\Seurat\third\_party\\ion>ion\build.bat -c opt_x64
```
2. Edit Ion's ZLib configuration header to disable <unistd.h>. Locate
   `third_party/ion/third_party/zlib/src/zconf.h`, ~line 434, and disable the
   definition of `Z_HAVE_UNISTD_H`.

3. Build the pipeline executable with:
```
c:\Seurat>bazel --bazelrc tools\bazel_windows.rc build -c opt seurat/pipeline:seurat
```
   Note that "opt" builds can be 100x faster and are necessary for the
   pipeline to complete in reasonable time.

4. Build the viewer executable with:
```
c:\Seurat>bazel --bazelrc tools\bazel_windows.rc build -c opt seurat/viewer/butterfly
```

### Supported Toolchain

#### Bazel 0.7.0

The build is verified to work with Bazel release 0.7.0. Other versions,
including newer ones, might not work.

#### MSVC 2015 Update 3

The build is verified to work with MSVC 2015 Update 3. Other versions,
including newer ones, might not work.
