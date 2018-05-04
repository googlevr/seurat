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

#ifndef VR_SEURAT_BASE_PROGRESS_H_
#define VR_SEURAT_BASE_PROGRESS_H_

#include <iostream>
#include <memory>
#include <mutex>  // NOLINT
#include <string>

namespace seurat {
namespace base {

// A progress reporter. Call GetProgress() to get the singleton instance.
// Progress reporting is disabled by default. Call GetProgress()->Enable(true)
// to turn it on.
//
// Progress reporter is safe to use in parallel_for loops, where multiple
// workers increment a range cooperatively. It is NOT safe to use by multiple
// threads that each want to report their own ranges.
//
// Progress can report single "events" or "ranges". For events, the description
// is output to the stream once.

// For ranges, a progress bar is output. It gets printed on BeginRange() and
// EndRange(). IncrementRange() only re-prints the progress bar if it has
// visually changed. It is okay to call EndRange() before the range is complete.
// The progress bar will be left partially filled in this case and terminated by
// a newline.
class Progress {
 public:
  class Timer {
   public:
    virtual ~Timer() = default;
    virtual void Reset() = 0;
    virtual double GetInS() const = 0;
  };

  Progress();
  Progress(int terminal_width, std::ostream* stream,
           std::shared_ptr<Timer> timer);

  // Turns progress reporting on or off.
  void Enable(bool enabled);

  // Report a single progress event. Progress must not be in a range when
  // this function is called. Otherwise it check fails.
  void Event(const std::string& description);

  // Start reporting a new range, i.e. an event with multiple steps.
  void BeginRange(const std::string& description, int total_steps);

  // Increments the current range by |steps|. Check fails if no range is open.
  void IncrementRange(int steps);

  // Ends the current range. Check fails if no range is open.
  void EndRange();

 private:
  // Prints the progress bar.
  void PrintProgressBar();

  // Resets range status to be not active.
  void ResetRange();

  // Whether progress reporting is enabled.
  bool enabled_;

  // The stream used for reporting progress.
  std::ostream& stream_;

  // Width of the terminal in characters. Used to determine progress bar size.
  int terminal_width_;

  // Mutex for accessing the output stream.
  std::mutex mutex_;

  // A flag that indicates if we are currently in a range.
  bool range_active_;

  // Description of the current range.
  std::string range_description_;

  // Total steps in the current range.
  int range_total_steps_;

  // Number of steps completed in the current range.
  int range_completed_steps_;

  // Number of slots in the progress bar that are filled. This is used to
  // determine if the progress bar needs a redraw.
  int range_slots_filled_;

  // Measuring the elapsed time since the last BeginRange.
  std::shared_ptr<Timer> timer_;
};

// A utility class for scoped progress ranges. The destructor calls EndRange
// automatically, making this helper safe for reporting in tasks with early
// return statements.
class ScopedProgressRange {
 public:
  ScopedProgressRange(const std::string& description, int total_steps);
  ~ScopedProgressRange();
  void IncrementRange(int increment);

 private:
  Progress* progress_;
};

// Get the singleton instance.
Progress* GetProgress();

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_PROGRESS_H_
