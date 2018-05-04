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

"""Defines for the opensource build of seurat."""

def compiled_zipassets(name, asset_file, srcs, deps=[], visibility=None):
  """Generates a C++ source file defining assets for a library or program.
  """
  compiled_file = asset_file.replace(".iad", ".cc")
  asset_path = "$(location " + asset_file + ")"
  if asset_file not in srcs:
    srcs += [asset_file]
  native.genrule(
      name = "%s_compiled" % name,
      srcs = srcs,
      outs = [compiled_file],
      cmd = "$(location @ion//:zipasset_generator) {0} $$(dirname {0}) $(OUTS)".format(asset_path),
      tools = ["@ion//:zipasset_generator"],
      message = "Generating zipasset file %s" % asset_file,
      visibility = visibility,
  )
  native.cc_library(
      name = name,
      srcs = [compiled_file],
      deps = [
          "@ion//:ionbase",
      ] + deps,
      visibility = visibility,
  )

def cc_test(name, linkopts = [], tags = [], **kwargs):
  """Constructs the test rules for a set of unittests
  """
  native.cc_test(
      name = name,
      linkopts = linkopts,
      linkstatic = 1,
      tags = tags,
      **kwargs
  )

def test_suite(name, tests = [], tags = []):
  """Constructs a test suite for a list of unittests.
  """
  native.test_suite(
      name = name,
      tests = tests,
      tags = tags,
  )

def portable_proto_library():
  """Not used.
  """
  pass

def init_ion_settings():
  """Not used.
  """
  pass

def select_for_ion_platform():
  """Not used.
  """
  pass
