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

#include "seurat/base/progress.h"

#include <cstdio>
#include <iostream>
#if defined(ION_PLATFORM_LINUX)
#include <sys/ioctl.h>
#include <unistd.h>
#elif defined(ION_PLATFORM_WINDOWS)
#include <windows.h>
#endif

#include "ion/base/logging.h"
#include "ion/base/staticsafedeclare.h"
#include "ion/port/timer.h"

// Removes Ion's snprintf macro, if present, since it interferes with use of
// std::snprintf.
#ifdef snprintf
#undef snprintf
#endif

namespace seurat {
namespace base {

namespace {

// Returns the width in characters of the stdout terminal.
#if defined(ION_PLATFORM_LINUX)
int GetTerminalWidth() {
  struct winsize size;
  if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &size) < 0) {
    return 80;
  }
  return size.ws_col;
}
#elif defined(ION_PLATFORM_WINDOWS)
int GetTerminalWidth() {
  HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
  if (handle == INVALID_HANDLE_VALUE) {
    return 80;
  }
  CONSOLE_SCREEN_BUFFER_INFO console_screen_buffer_info;
  GetConsoleScreenBufferInfo(handle, &console_screen_buffer_info);
  return console_screen_buffer_info.dwSize.X;
}
#else
int GetTerminalWidth() { return 80; }
#endif

class IonTimer : public Progress::Timer {
 public:
  void Reset() override { ion_timer_.Reset(); }
  double GetInS() const override { return ion_timer_.GetInS(); }

 private:
  ion::port::Timer ion_timer_;
};

}  // namespace

Progress::Progress()
    : Progress(GetTerminalWidth(), &std::cout, std::make_shared<IonTimer>()) {}

Progress::Progress(int terminal_width, std::ostream* stream,
                   std::shared_ptr<Timer> timer)
    : enabled_(false),
      stream_(*stream),
      terminal_width_(terminal_width),
      timer_(std::move(timer)) {
  ResetRange();
}

void Progress::Enable(bool enabled) {
  std::lock_guard<std::mutex> lock(mutex_);
  enabled_ = enabled;
}

void Progress::ResetRange() {
  range_active_ = false;
  range_description_ = std::string();
  range_total_steps_ = 0;
  range_completed_steps_ = 0;
  // Set to -1 so that we also draw the empty range on the first call to
  // PrintProgressBar.
  range_slots_filled_ = -1;
}

void Progress::Event(const std::string& description) {
  std::lock_guard<std::mutex> lock(mutex_);
  if (!enabled_) {
    return;
  }
  CHECK(!range_active_);
  stream_ << description << "\n";
}

void Progress::BeginRange(const std::string& description, int total_steps) {
  std::lock_guard<std::mutex> lock(mutex_);
  if (!enabled_) {
    return;
  }
  CHECK(!range_active_);
  range_active_ = true;
  range_description_ = description;
  range_total_steps_ = total_steps;
  range_completed_steps_ = 0;
  range_slots_filled_ = -1;
  timer_->Reset();
  PrintProgressBar();
}

void Progress::IncrementRange(int steps) {
  std::lock_guard<std::mutex> lock(mutex_);
  if (!enabled_) {
    return;
  }
  CHECK(range_active_);
  range_completed_steps_ += steps;
  CHECK_LE(range_completed_steps_, range_total_steps_);
  PrintProgressBar();
}

void Progress::EndRange() {
  std::lock_guard<std::mutex> lock(mutex_);
  if (!enabled_) {
    return;
  }
  CHECK(range_active_);
  PrintProgressBar();
  stream_ << "\n";
  ResetRange();
}

void Progress::PrintProgressBar() {
  int raw_seconds = static_cast<int>(timer_->GetInS());
  int hours = (raw_seconds / 3600) % 24;
  int minutes = (raw_seconds / 60) % 60;
  int seconds = raw_seconds % 60;
  char buffer[9];
  std::snprintf(buffer, sizeof(buffer), "%02d:%02d:%02d", hours, minutes,
                seconds);
  std::string elapsed_time = buffer;

  // -5 for ": [] ", and -1 on Windows for some unknown reason.
  const int bar_size =
      terminal_width_ - range_description_.size() - elapsed_time.size() - 5 - 1;

  // If we don't have enough space for the bar, then print the description only.
  if (bar_size <= 0) {
    stream_ << range_description_;
    return;
  }

  const int slots_filled =
      range_completed_steps_ * bar_size / range_total_steps_;
  // Don't redraw, if no new slots were filled.
  if (slots_filled == range_slots_filled_) {
    return;
  }
  range_slots_filled_ = slots_filled;
  const int slots_empty = bar_size - slots_filled;

  stream_ << range_description_ << ": [" << std::string(slots_filled, '+')
          << std::string(slots_empty, ' ') << "] " << elapsed_time << "\r";
  stream_.flush();
}

ScopedProgressRange::ScopedProgressRange(const std::string& description,
                                         int total_steps)
    : progress_(GetProgress()) {
  progress_->BeginRange(description, total_steps);
}

ScopedProgressRange::~ScopedProgressRange() { progress_->EndRange(); }

void ScopedProgressRange::IncrementRange(int increment) {
  progress_->IncrementRange(increment);
}

Progress* GetProgress() {
  ION_DECLARE_SAFE_STATIC_POINTER(Progress, progress);
  return progress;
}

}  // namespace base
}  // namespace seurat
