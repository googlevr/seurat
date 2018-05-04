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

#include <condition_variable>  // NOLINT(build/c++11)
#include <memory>
#include <mutex>  // NOLINT(build/c++11)
#include <vector>

#include "ion/base/logging.h"
#include "ion/base/staticsafedeclare.h"
#include "absl/synchronization/barrier.h"
#include "seurat/base/util.h"

namespace seurat {
namespace base {
namespace {

// This class encapsulates a thread that can run tasks passed to it.  The thread
// is started on creation, and is joined on destruction of the class instance.
class ThreadedRunner {
 public:
  ThreadedRunner()
      : is_running_(true),
        func_(nullptr),
        barrier_(nullptr),
        count_(0),
        thread_(&ThreadedRunner::ThreadFunc, this) {}

  ~ThreadedRunner() {
    {
      std::unique_lock<std::mutex> lock(mutex_);
      is_running_ = false;
      condvar_.notify_all();
    }
    thread_.join();
  }

  // Run |func| on the encapsulated thread.  Returns true iff the |func| was
  // successfully scheduled; may fail if this ThreadedRunner is already busy.
  // The task is run for an index starting at |thread_index|, up to |count|,
  // incrementing by |thread_count|.
  bool Run(std::function<void(int)>* func, absl::Barrier* barrier,
           int thread_index, int thread_count, int count) {
    std::unique_lock<std::mutex> lock(mutex_);
    if (func_ != nullptr) {
      return false;
    }
    func_ = func;
    barrier_ = barrier;
    thread_index_ = thread_index;
    thread_count_ = thread_count;
    count_ = count;
    condvar_.notify_all();
    return true;
  }

 private:
  // Thread main function.  This function loops waiting for tasks and executing
  // them, until it is requested to quit by setting |is_running_| to false.
  void ThreadFunc() {
    std::unique_lock<std::mutex> lock(mutex_);
    while (true) {
      // Wait for an incoming task.
      while (is_running_ && func_ == nullptr) {
        condvar_.wait(lock);
      }
      if (!is_running_) {
        return;
      }
      DCHECK(func_ != nullptr);

      // Execute!
      lock.unlock();
      for (int i = thread_index_; i < count_; i += thread_count_) {
        (*func_)(i);
      }
      if (barrier_->Block()) delete barrier_;

      // Mark this ThreadedRunner as available for another task.
      lock.lock();
      func_ = nullptr;
      barrier_ = nullptr;
      thread_index_ = 0;
      thread_count_ = 0;
      count_ = 0;
    }
  }

  // ThreadedRunner execution state.

  // Protects state shared between Run() in the public interface and
  // ThreadFunc() executing the thread logic.
  std::mutex mutex_;

  // Notified when execution state has changed (e.g. a new task is added, or the
  // thread is asked to quit.
  std::condition_variable condvar_;

  // True iff the thread should continue running.
  bool is_running_;

  // ThreadedRunner per-task state.  Set every time a task is scheduled, and
  // reset when it is completed.

  // The task to run.  nullptr iff the thread is idle.
  std::function<void(int)>* func_;

  // The Barrier to wait on task completion.  nullptr iff the thread is idle.
  absl::Barrier* barrier_;

  // The index of this thread in running |func_|.  Used when distributing
  // execution of |func_| across many ThreadedRunner instances.
  int thread_index_;

  // The total number of threads in this run of |func_|.  Used to increment the
  // index for each ThreadedRunner, when distributing execution of |func_|
  // across many ThreadedRunner instances.
  int thread_count_;

  // The total count of executions of |func_| requested, across all
  // ThreadedRunner instances.
  int count_;

  // The std::thread is last, since it has to be initialized after all other
  // members are initialized.
  std::thread thread_;
};

// This class maintains a pool of ThreadedRunner instances to run tasks.  More
// instances are allocated as necessary, and reused if possible.
class ThreadedRunnerPool {
 public:
  ThreadedRunnerPool() {}

  // Run |func| using this ThreadedRunnerPool.
  void Run(int thread_count, int count, std::function<void(int)>* func) {
    absl::Barrier* barrier = new absl::Barrier(thread_count);

    // Run |thread_count| - 1 tasks on worker threads.
    {
      std::unique_lock<std::mutex> lock(mutex_);
      int started = 0;
      int index = 0;
      while (started + 1 < thread_count) {
        if (index >= runners_.size()) {
          runners_.emplace_back(new ThreadedRunner());
        }
        if (runners_[index]->Run(func, barrier, started + 1, thread_count,
                                 count)) {
          ++started;
        }
        ++index;
      }
    }

    // Run one task on the main thread.
    for (int i = 0; i < count; i += thread_count) {
      (*func)(i);
    }
    if (barrier->Block()) delete barrier;
  }

 private:
  // Protects access to |runners_|.
  std::mutex mutex_;

  // The ThreadedRunners running tasks in this ThreadedRunnerPool.
  std::vector<std::unique_ptr<ThreadedRunner>> runners_;
};

}  // namespace

void ParallelFor(int thread_count, int count, std::function<void(int)> func) {
  DCHECK_GT(thread_count, 0);
  ION_DECLARE_SAFE_STATIC_POINTER(ThreadedRunnerPool, s_Pool);
  s_Pool->Run(thread_count, count, &func);
}

void BalancedParallelFor(int thread_count, int count,
                         std::function<void(int)> func) {
  std::atomic<int> counter(0);
  ParallelFor(thread_count, thread_count, [&](int tid) {
    while (true) {
      int index = counter.fetch_add(1);
      if (index >= count) return;
      func(index);
    }
  });
}

}  // namespace base
}  // namespace seurat
