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

#ifndef VR_SEURAT_BASE_BYTESTREAM_H_
#define VR_SEURAT_BASE_BYTESTREAM_H_

#include <cstddef>
#include <string>

#include "absl/strings/string_view.h"

namespace seurat {
namespace base {

// A ByteSink is an abstract interface for an object that consumes a sequence of
// bytes.
class ByteSink {
 public:
  ByteSink() {}
  ByteSink(const ByteSink&) = delete;
  ByteSink& operator=(const ByteSink&) = delete;
  virtual ~ByteSink() {}

  // Appends "n" bytes from "bytes" into the ByteSink.
  virtual void Append(const char* bytes, size_t n) = 0;
};

class StringByteSink : public ByteSink {
 public:
  explicit StringByteSink(std::string* dest) : dest_(dest) {}
  StringByteSink(const StringByteSink&) = delete;
  StringByteSink& operator=(const StringByteSink&) = delete;

  void Append(const char* data, size_t n) override { dest_->append(data, n); }

 private:
  std::string* dest_;
};

class UncheckedArrayByteSink : public ByteSink {
 public:
  explicit UncheckedArrayByteSink(char* dest) : dest_(dest) {}
  UncheckedArrayByteSink(const UncheckedArrayByteSink&) = delete;
  UncheckedArrayByteSink& operator=(const UncheckedArrayByteSink&) = delete;
  void Append(const char* data, size_t n) override;

 private:
  char* dest_;
};

// A ByteSource is an abstract interface for an object that produces a sequence
// of bytes.
class ByteSource {
 public:
  ByteSource() {}
  ByteSource(const ByteSource&) = delete;
  ByteSource& operator=(const ByteSource&) = delete;
  virtual ~ByteSource() {}

  // Returns the number of bytes left to read from the source. Available()
  // should decrease by N each time Skip(N) is called; Available() may not
  // increase. If Available() returns 0, that indicates that the ByteSource is
  // exhausted.
  virtual size_t Available() const = 0;

  // Returns an absl::string_view of the next contiguous region of the source.
  // Does not reposition the source. The returned region is empty iff
  // Available() == 0.
  //
  // The returned region is valid until the next call to Skip() or until this
  // object is destroyed, whichever occurs first.
  //
  // The length of the returned absl::string_view will be <= Available().
  virtual absl::string_view Peek() = 0;

  // Skips the next n bytes. Invalidates any absl::string_view returned by a
  // previous call to Peek().
  //
  // REQUIRES: Available() >= n
  virtual void Skip(size_t n) = 0;

  // Writes the next n bytes in this ByteSource to the given ByteSink, and
  // advances this ByteSource past the copied bytes. The default implementation
  // of this method just copies the bytes normally, but subclasses might
  // override CopyTo to optimize certain cases.
  //
  // REQUIRES: Available() >= n
  virtual void CopyTo(ByteSink* sink, size_t n);
};

class ArrayByteSource : public ByteSource {
 public:
  explicit ArrayByteSource(absl::string_view s) : input_(s) {}
  ArrayByteSource(const ArrayByteSource&) = delete;
  ArrayByteSource& operator=(const ArrayByteSource&) = delete;

  size_t Available() const override { return input_.size(); }

  absl::string_view Peek() override { return input_; }

  void Skip(size_t n) override {
    assert(n <= input_.size());
    input_.remove_prefix(n);
  }

 private:
  absl::string_view input_;
};

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_BYTESTREAM_H_
