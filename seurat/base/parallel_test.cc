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

#include "seurat/base/parallel.h"

#include <array>
#include <atomic>
#include <condition_variable>  // NOLINT(build/c++11)
#include <mutex>               // NOLINT(build/c++11)

#include "gtest/gtest.h"

namespace seurat {
namespace base {
namespace {

// Tests multi-threaded increment of a single variable, using increasing numbers
// of threads.
TEST(ParallelTest, IncreasingThreadCounts) {
  for (int thread_count = 1; thread_count < 16; ++thread_count) {
    std::mutex mutex;
    std::condition_variable condvar;
    int index = 0;
    ParallelFor(thread_count, 1000, [&](int i) {
      std::unique_lock<std::mutex> lock(mutex);
      while (i != index) {
        condvar.wait(lock);
      }
      ++index;
      condvar.notify_all();
    });
    EXPECT_EQ(1000, index);
  }
}

// Tests recursive use of ParallelFor().
TEST(ParallelTest, RecursiveParallelFor) {
  std::atomic<int> count(0);
  ParallelFor(16, 16, [&](int i) {
    ParallelFor(16, 16, [&](int j) {
      for (int k = 0; k < 16; ++k) {
        ++count;
      }
    });
  });
  EXPECT_EQ(16 * 16 * 16, count.load());
}

TEST(ParallelTest, BalancedParallelFor) {
  const int kNumElements = 16 * 4 + 7;
  std::array<std::atomic<int>, kNumElements> counters;
  for (auto& c : counters) {
    c = 0;
  }
  BalancedParallelFor(16, kNumElements, [&](int i) { counters[i]++; });
  for (const auto& count : counters) {
    EXPECT_EQ(1, count);
  }
}

}  // namespace
}  // namespace base
}  // namespace seurat
