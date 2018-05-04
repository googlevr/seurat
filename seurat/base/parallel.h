/*
Copyright 2017 Google Inc. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS-IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef VR_SEURAT_BASE_PARALLEL_H_
#define VR_SEURAT_BASE_PARALLEL_H_

#include <algorithm>
#include <functional>
#include <thread>  // NOLINT(build/c++11)

namespace seurat {
namespace base {

// Returns the number of hardware threads, always returns more than zero.
inline int GetNumberOfHardwareThreads() {
  // According to the standard, hardware_concurrency() may return zero if "this
  // value is not computable or well defined". In that case, we want to have a
  // thread count of one (instead of zero).
  return std::max(static_cast<unsigned int>(1),
                  std::thread::hardware_concurrency());
}

// Repeatedly executes |func| multiple times, specified by |count|.
// Different executions of |func| may occur on different threads.
//
// Func has function signature void(int i), with |i| taking values in the range
// [0, |count|).
void ParallelFor(int thread_count, int count, std::function<void(int)> func);

// Like ParallelFor but automatically balances execution.
//
// Instead of executing ~count/thread_count consecutive invocations of |func| on
// the same thread, each thread will grab the next index from a shared counter.
void BalancedParallelFor(int thread_count, int count,
                         std::function<void(int)> func);

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_PARALLEL_H_
