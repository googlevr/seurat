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
#include <sstream>

#include "gtest/gtest.h"

namespace seurat {
namespace base {
namespace {

class FakeTimer : public Progress::Timer {
 public:
  FakeTimer() { seconds_ = 0.0; }
  void Reset() override { seconds_ = 0.0; }
  double GetInS() const override { return seconds_; }
  void Set(double seconds) { seconds_ = seconds; }

 private:
  double seconds_;
};

TEST(Progress, Disabled) {
  const std::string kEvent("Test");
  std::ostringstream stream;
  const int kTerminalWidth = 80;
  Progress progress(kTerminalWidth, &stream,
                    std::unique_ptr<Progress::Timer>(new FakeTimer()));
  // Progress reporting is disabled by default. Test that it doesn't output any
  // progress reports.
  progress.Event(kEvent);
  EXPECT_EQ(std::string(), stream.str());
  progress.BeginRange("TestRange", 1);
  EXPECT_EQ(std::string(), stream.str());
  progress.IncrementRange(1);
  EXPECT_EQ(std::string(), stream.str());
  progress.EndRange();
  EXPECT_EQ(std::string(), stream.str());
}

TEST(Progress, BeginRangeInsideRangeChecks) {
  std::ostringstream stream;
  const int kTerminalWidth = 80;
  Progress progress(kTerminalWidth, &stream,
                    std::unique_ptr<Progress::Timer>(new FakeTimer()));
  progress.Enable(true);
  progress.BeginRange("Test", 1);
  EXPECT_DEATH(progress.BeginRange("Test", 1), "CHECK");
}

TEST(Progress, EventInsideRangeChecks) {
  std::ostringstream stream;
  const int kTerminalWidth = 80;
  Progress progress(kTerminalWidth, &stream,
                    std::unique_ptr<Progress::Timer>(new FakeTimer()));
  progress.Enable(true);
  progress.BeginRange("Test range", 1);
  EXPECT_DEATH(progress.Event("Test event"), "CHECK");
}

TEST(Progress, IncrementOutsideRangeChecks) {
  std::ostringstream stream;
  const int kTerminalWidth = 80;
  Progress progress(kTerminalWidth, &stream,
                    std::unique_ptr<Progress::Timer>(new FakeTimer()));
  progress.Enable(true);
  EXPECT_DEATH(progress.IncrementRange(1), "CHECK");
}

TEST(Progress, IncrementBeyonEndChecks) {
  std::ostringstream stream;
  const int kTerminalWidth = 80;
  Progress progress(kTerminalWidth, &stream,
                    std::unique_ptr<Progress::Timer>(new FakeTimer()));
  progress.Enable(true);
  progress.BeginRange("Test", 1);
  progress.IncrementRange(1);
  EXPECT_DEATH(progress.IncrementRange(1), "CHECK");
}

TEST(Progress, EndRangeOutsideRangeChecks) {
  std::ostringstream stream;
  const int kTerminalWidth = 80;
  Progress progress(kTerminalWidth, &stream,
                    std::unique_ptr<Progress::Timer>(new FakeTimer()));
  progress.Enable(true);
  EXPECT_DEATH(progress.EndRange(), "CHECK");
}

TEST(Progress, EndRangeBeforeEndWorks) {
  std::ostringstream stream;
  const int kTerminalWidth = 80;
  Progress progress(kTerminalWidth, &stream,
                    std::unique_ptr<Progress::Timer>(new FakeTimer()));
  progress.Enable(true);
  progress.BeginRange("Test", 1);
  stream.str(std::string());
  progress.EndRange();
  EXPECT_EQ("\n", stream.str());
}

TEST(Progress, Event) {
  const std::string kEvent("Test");
  std::ostringstream stream;
  const int kTerminalWidth = 80;
  Progress progress(kTerminalWidth, &stream,
                    std::unique_ptr<Progress::Timer>(new FakeTimer()));
  progress.Enable(true);
  progress.Event(kEvent);
  EXPECT_EQ(kEvent + "\n", stream.str());
}

TEST(Progress, Range) {
  const std::string kEvent("Test");
  std::ostringstream stream;
  // 10 for the steps, 4 for "Test", 4 for ": []" and 1 for Windows.
  // This width results in exactly one "+" per step.
  const int kTerminalWidth = 28;
  const int kSteps = 10;
  std::shared_ptr<FakeTimer> timer = std::make_shared<FakeTimer>();
  Progress progress(kTerminalWidth, &stream, timer);
  progress.Enable(true);

  progress.BeginRange(kEvent, kSteps);
  std::string expected_bar(kEvent + ": [          ] 00:00:00\r");
  EXPECT_EQ(expected_bar, stream.str());
  // Clear the stream after begin and every step. Otherwise it accumulates all
  // stages of the progress bar.
  stream.str(std::string());

  for (int i = 1; i <= kSteps; ++i) {
    timer->Set(i);
    progress.IncrementRange(1);
    char expected_elapsed_time[9];
    std::snprintf(expected_elapsed_time, sizeof(expected_elapsed_time),
                  "00:00:%02d", i);
    expected_bar = std::string(kEvent + ": [" + std::string(i, '+') +
                               std::string(kSteps - i, ' ') + "] " +
                               expected_elapsed_time + "\r");
    EXPECT_EQ(expected_bar, stream.str());
    stream.str(std::string());
  }

  progress.EndRange();
  expected_bar = "\n";
  EXPECT_EQ(expected_bar, stream.str());
}

TEST(Progress, RangeWithLongDescription) {
  const std::string kEvent("TestWithLongDesciption");
  std::ostringstream stream;
  // Not wide enough for kEvent and progress bar. Expect that only kEvent is
  // printed.
  const int kTerminalWidth = 10;
  const int kSteps = 10;
  Progress progress(kTerminalWidth, &stream,
                    std::unique_ptr<Progress::Timer>(new FakeTimer()));
  progress.Enable(true);

  progress.BeginRange(kEvent, kSteps);
  EXPECT_EQ(kEvent, stream.str());
  // Clear the stream after begin and every step. Otherwise it accumulates all
  // stages of the progress bar.
  stream.str(std::string());

  for (int i = 1; i <= kSteps; ++i) {
    progress.IncrementRange(1);
    EXPECT_EQ(kEvent, stream.str());
    stream.str(std::string());
  }

  progress.EndRange();
  EXPECT_EQ(std::string("TestWithLongDesciption\n"), stream.str());
}

}  // namespace
}  // namespace base
}  // namespace seurat
