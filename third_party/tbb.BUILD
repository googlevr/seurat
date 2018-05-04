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
    name = "windows",
    values = {
        "crosstool_top": "//crosstools/windows",
    },
)

cc_library(
    name = "tbb",
    srcs = [
        "include/tbb/aligned_space.h",
        "include/tbb/atomic.h",
        "include/tbb/cache_aligned_allocator.h",
        "include/tbb/concurrent_hash_map.h",
        "include/tbb/concurrent_vector.h",
        "include/tbb/critical_section.h",
        "include/tbb/internal/_concurrent_queue_impl.h",
        "include/tbb/internal/_concurrent_unordered_impl.h",
        "include/tbb/machine/gcc_ia32_common.h",
        "include/tbb/machine/gcc_itsx.h",
        "include/tbb/machine/linux_common.h",
        "include/tbb/machine/linux_intel64.h",
        "include/tbb/partitioner.h",
        "include/tbb/pipeline.h",
        "include/tbb/queuing_mutex.h",
        "include/tbb/queuing_rw_mutex.h",
        "include/tbb/reader_writer_lock.h",
        "include/tbb/recursive_mutex.h",
        "include/tbb/spin_mutex.h",
        "include/tbb/spin_rw_mutex.h",
        "include/tbb/task_scheduler_observer.h",
        "include/tbb/tbb_allocator.h",
        "include/tbb/tbb_exception.h",
        "include/tbb/tbb_machine.h",
        "include/tbb/tbb_profiling.h",
        "include/tbb/tbb_stddef.h",
        "include/tbb/tick_count.h",
        "src/old/concurrent_vector_v2.cpp",
        "src/old/concurrent_vector_v2.h",
        "src/old/spin_rw_mutex_v2.cpp",
        "src/old/spin_rw_mutex_v2.h",
        "src/old/task_v2.cpp",
        "src/rml/client/rml_factory.h",
        "src/rml/client/rml_tbb.cpp",
        "src/rml/include/rml_base.h",
        "src/rml/include/rml_tbb.h",
        "src/rml/server/thread_monitor.h",
        "src/tbb/arena.cpp",
        "src/tbb/arena.h",
        "src/tbb/cache_aligned_allocator.cpp",
        "src/tbb/cilk-tbb-interop.h",
        "src/tbb/concurrent_hash_map.cpp",
        "src/tbb/concurrent_monitor.cpp",
        "src/tbb/concurrent_monitor.h",
        "src/tbb/concurrent_queue.cpp",
        "src/tbb/concurrent_vector.cpp",
        "src/tbb/condition_variable.cpp",
        "src/tbb/critical_section.cpp",
        "src/tbb/custom_scheduler.h",
        "src/tbb/dynamic_link.cpp",
        "src/tbb/dynamic_link.h",
        "src/tbb/governor.cpp",
        "src/tbb/governor.h",
        "src/tbb/intrusive_list.h",
        "src/tbb/itt_notify.cpp",
        "src/tbb/itt_notify.h",
        "src/tbb/mailbox.h",
        "src/tbb/market.cpp",
        "src/tbb/market.h",
        "src/tbb/mutex.cpp",
        "src/tbb/observer_proxy.cpp",
        "src/tbb/observer_proxy.h",
        "src/tbb/pipeline.cpp",
        "src/tbb/private_server.cpp",
        "src/tbb/queuing_mutex.cpp",
        "src/tbb/queuing_rw_mutex.cpp",
        "src/tbb/reader_writer_lock.cpp",
        "src/tbb/recursive_mutex.cpp",
        "src/tbb/scheduler.cpp",
        "src/tbb/scheduler.h",
        "src/tbb/scheduler_common.h",
        "src/tbb/scheduler_utility.h",
        "src/tbb/semaphore.cpp",
        "src/tbb/semaphore.h",
        "src/tbb/spin_mutex.cpp",
        "src/tbb/spin_rw_mutex.cpp",
        "src/tbb/task.cpp",
        "src/tbb/task_group_context.cpp",
        "src/tbb/task_stream.h",
        "src/tbb/tbb_assert_impl.h",
        "src/tbb/tbb_main.cpp",
        "src/tbb/tbb_main.h",
        "src/tbb/tbb_misc.cpp",
        "src/tbb/tbb_misc.h",
        "src/tbb/tbb_misc_ex.cpp",
        "src/tbb/tbb_statistics.cpp",
        "src/tbb/tbb_statistics.h",
        "src/tbb/tbb_thread.cpp",
        "src/tbb/tbb_version.h",
        "src/tbb/tls.h",
    ],
    hdrs = glob([
        "include/tbb/blocked_range.h",
        "include/tbb/compat/condition_variable",
        "include/tbb/concurrent_unordered_map.h",
        "include/tbb/concurrent_unordered_set.h",
        "include/tbb/mutex.h",
        "include/tbb/parallel_for.h",
        "include/tbb/parallel_reduce.h",
        "include/tbb/task.h",
        "include/tbb/task_scheduler_init.h",
        "include/tbb/tbb_thread.h",
        "include/tbb/*.h",
        "include/tbb/internal/*.h",
    ]),
    copts = select({
        ":windows": ["-DUSE_WINTHREAD"],
        "//conditions:default": ["-DUSE_PTHREAD"],
    }) + [
        # "-DDO_ITT_NOTIFY",    # We don't use intels ITT.
        "-DTBB_PREVIEW_GRAPH_NODES=1",
        "-D__TBB_BUILD=1",
        "-w",  # Inhibits warnings.
        "-fPIC",
        "-Wno-parentheses",
        "-Wno-non-virtual_dtor",
    ],
    includes = [
        "include",
        "src",
        "src/rml/include",
        "tbb/include",
    ],
    linkopts =
        select({
            ":windows": [],
            "//conditions:default": [
                "-ldl",
                "-lpthread",
                "-lrt",
            ],
        }),
    deps = [":tbb_ver"],
)

# Generate version_string.ver which is required by tbb_version.h
cc_library(
    name = "tbb_ver",
    srcs = [],
    hdrs = ["version_string.ver"],
)

VERSION_STRING_VER = """
#define __TBB_VERSION_STRINGS(N) "Empty"
"""

genrule(
    name = "generate_version_string",
    outs = ["version_string.ver"],
    cmd = ("echo '%s' > $(OUTS)" % VERSION_STRING_VER),
    visibility = ["//visibility:private"],
)
