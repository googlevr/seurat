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

#include "seurat/base/structured_io.h"

#include <memory>
#include <random>
#include <string>
#include <vector>

#include "ion/gfx/image.h"
#include "ion/math/matrix.h"
#include "ion/math/vector.h"
#include "gtest/gtest.h"
#include "absl/strings/string_view.h"
#include "seurat/base/array2d.h"
#include "seurat/base/array2d_util.h"
#include "seurat/base/array2d_view.h"
#include "seurat/base/color.h"
#include "seurat/base/ion_util_no_gl.h"
#include "seurat/base/structured_io_tests_utils.h"

namespace seurat {
namespace base {
namespace {

using ion::math::Matrix4f;
using ion::math::Point2i;
using ion::math::Vector2i;

// The number of test values (of random types) written to the ByteSink for
// testing.
static const int kFuzzTestValueCount = 32 * 1024;
using PRNG = std::mt19937;

// Tests a simple stream of POD values.
TEST(StructuredIoTest, PodValuesTest) {
  std::string buf;
  base::StringByteSink byte_sink(&buf);
  StructureSink sink(&byte_sink);
  sink.WritePod<int>(-1);
  sink.WritePod<unsigned char>(2);
  sink.WritePod<uint64>(~0);
  sink.WritePod<double>(13.37);
  const int sized_array[3] = {1, 2, 3};
  sink.WritePodArray(sized_array);

  const std::vector<char> message{'H', 'e', 'l', 'l', 'o', '!', '\0'};
  const int kDataSize = 7;
  CHECK(kDataSize == message.size())
      << "Constant required for array readback variable declaration "
         "(output_unsized_array), and constant must match source data "
         "array size (message).";
  sink.WritePodArray(message.data(), kDataSize);

  const int sized_3d_array[2][3][2] = {{{1, 2}, {3, 4}, {8, 14}},
                                       {{5, 6}, {314, 100}, {13, 11}}};
  sink.WritePodArray(sized_3d_array);

  ArrayByteSource byte_source(absl::string_view(buf.data(), buf.size()));
  StructureSource source(&byte_source);
  EXPECT_EQ(-1, source.ReadPod<int>());
  EXPECT_EQ(2, source.ReadPod<unsigned char>());
  EXPECT_EQ(static_cast<uint64>(~0), source.ReadPod<uint64>());
  EXPECT_EQ(13.37, source.ReadPod<double>());

  int output_array[3] = {0};
  source.ReadPodArray(output_array);
  EXPECT_EQ(sized_array[0], output_array[0]);
  EXPECT_EQ(sized_array[1], output_array[1]);
  EXPECT_EQ(sized_array[2], output_array[2]);

  char output_unsized_array[kDataSize] = {0};
  source.ReadPodArray(output_unsized_array, kDataSize);
  const int indices[] = {0, 1, 2, 3, 4, 5, 6};
  static_assert(ARRAYSIZE(indices) == kDataSize, "Please test all indices.");
  for (int i : indices) {
    EXPECT_EQ(output_unsized_array[i], message[i]);
  }

  int output_3d_array[2][3][2] = {{{0}}};
  static_assert(sizeof(output_3d_array) == sizeof(sized_3d_array),
                "In/out 3D array sized must match.");
  source.ReadPodArray(output_3d_array);
  for (int outer_index = 0; outer_index < ARRAYSIZE(output_3d_array);
       ++outer_index) {
    for (int middle_index = 0; middle_index < ARRAYSIZE(output_3d_array[0]);
         ++middle_index) {
      for (int inner_index = 0; inner_index < ARRAYSIZE(output_3d_array[0][0]);
           ++inner_index) {
        EXPECT_EQ(sized_3d_array[outer_index][middle_index][inner_index],
                  output_3d_array[outer_index][middle_index][inner_index]);
      }
    }
  }
}

// Hides the NonPod type inside a namespace to ensure that
// argument-dependent-lookup works for WriteObjectsToSink and
// ReadObjectsFromSource.
namespace check_adl {

// Defines a class containing only plain-old-data, but violating POD constraints
// due to the constructor.
class NonPod {
 public:
  NonPod() : x_{} {}
  explicit NonPod(int x) : x_(x) {}

  bool operator==(const NonPod& rhs) const { return rhs.x_ == x_; }
  bool operator!=(const NonPod& rhs) const { return !operator==(rhs); }

  int x_;
};

static_assert(!std::is_pod<NonPod>::value,
              "NonPod must be non-POD to test ADL.");

inline void WriteObjectsToSink(const NonPod* objects, size_t count,
                               StructureSink* sink) {
  sink->WriteObjects(objects, count);
}

inline void ReadObjectsFromSource(StructureSource* source, NonPod* objects,
                                  size_t count) {
  source->ReadObjects(objects, count);
}

}  // namespace check_adl

// Tests structured IO adapters for non-POD values.
TEST(StructuredIoTest, NonPodAdapterTest) {
  using check_adl::NonPod;
  static_assert(
      !std::is_pod<NonPod>::value,
      "NonPodAdapterTest requires NonPod to not be a plain-old-data type.");

  std::string buf;
  base::StringByteSink byte_sink(&buf);
  StructureSink sink(&byte_sink);

  std::vector<NonPod> source_array{NonPod(1)};
  sink.WritePodArray(source_array);

  ArrayByteSource byte_source(absl::string_view(buf.data(), buf.size()));
  StructureSource source(&byte_source);

  std::vector<NonPod> final_array;
  source.ReadPodArray(&final_array);
  EXPECT_EQ(source_array, final_array);
  EXPECT_GT(final_array.size(), 0);
}

// Tests a simple stream of vector/point/color/matrix values.
TEST(StructuredIoTest, VectorMatrixValuesTest) {
  std::string buf;
  base::StringByteSink byte_sink(&buf);
  StructureSink sink(&byte_sink);

  const ion::math::Vector2f vector_value(512.66, 73.48);
  const ion::math::Point2i point_value(-255, 255);
  const Color4d color_value(0.1, 0.2, 0.3, 1.0);
  const ion::math::Matrix2d matrix_value(-1.0, -0.0, 0.0, 1.0);

  sink.WriteVector(vector_value);
  sink.WritePoint(point_value);
  sink.WriteColor(color_value);
  sink.WriteMatrix(matrix_value);

  ArrayByteSource byte_source(absl::string_view(buf.data(), buf.size()));
  StructureSource source(&byte_source);
  EXPECT_EQ(vector_value, source.ReadVector<ion::math::Vector2f>());
  EXPECT_EQ(point_value, source.ReadPoint<ion::math::Point2i>());
  EXPECT_EQ(color_value, source.ReadColor<Color4d>());
  EXPECT_EQ(matrix_value, source.ReadMatrix<ion::math::Matrix2d>());
}

// Tests a simple stream of string/buffer values.
TEST(StructuredIoTest, StringBufferValuesTest) {
  std::string buf;
  base::StringByteSink byte_sink(&buf);
  StructureSink sink(&byte_sink);

  const std::string string_value("This is a test string.");
  const std::vector<char> buffer_value{3, 1, 4, 1, 5, 2, 6, 5, 3, 5, 9};

  sink.WriteString(string_value);
  sink.WriteBytes(buffer_value.data(), buffer_value.size());

  ArrayByteSource byte_source(absl::string_view(buf.data(), buf.size()));
  StructureSource source(&byte_source);
  EXPECT_EQ(string_value, source.ReadString());
  std::vector<char> read_buffer(buffer_value.size());
  source.ReadBytes(read_buffer.data(), read_buffer.size());
  EXPECT_EQ(buffer_value, read_buffer);
}

TEST(StructuredIoTest, WriteAndReadImage) {
  const Vector2i kSize(32, 16);
  const ion::gfx::Image::Format kFormat = ion::gfx::Image::kRgba8888;

  // Build an Ion image from a Seurat Image.
  Array2D<Color4ui8> source_image(kSize);
  for (int y = 0; y < kSize[1]; ++y) {
    for (int x = 0; x < kSize[0]; ++x) {
      const Point2i p(x, y);
      source_image.At(p) = Color4ui8(x, y, x + y, 255);
    }
  }
  ion::gfx::ImagePtr ion_image = CreateImage(kFormat, source_image.GetSize());
  MutableArray2DView<Color4ui8> ion_image_view(
      ion_image->GetData()->GetMutableData<Color4ui8>(),
      Vector2i(ion_image->GetWidth(), ion_image->GetHeight()));
  CopyArray(source_image, &ion_image_view, Vector2i::Zero());

  // Serialize.
  std::string buf;
  base::StringByteSink byte_sink(&buf);
  StructureSink sink(&byte_sink);
  sink.WriteImage(ion_image);

  // Deserialize.
  ArrayByteSource byte_source(absl::string_view(buf.data(), buf.size()));
  StructureSource source(&byte_source);
  ion::gfx::ImagePtr deserialized_ion_image = source.ReadImage();

  // Validate.
  Array2DView<Color4ui8> deserialized_ion_image_view(
      deserialized_ion_image->GetData()->GetData<Color4ui8>(),
      Vector2i(deserialized_ion_image->GetWidth(),
               deserialized_ion_image->GetHeight()));

  EXPECT_EQ(kFormat, deserialized_ion_image->GetFormat());
  EXPECT_EQ(kSize[0], deserialized_ion_image->GetWidth());
  EXPECT_EQ(kSize[1], deserialized_ion_image->GetHeight());

  for (int y = 0; y < kSize[1]; ++y) {
    for (int x = 0; x < kSize[0]; ++x) {
      const Point2i p(x, y);
      EXPECT_EQ(deserialized_ion_image_view.At(p), source_image.At(p)) << p;
    }
  }
}

// Implements tag dispatch for IteratorFuncWithType.
template <typename T>
struct TypeParam {};

// Perform |iterations| iterations of picking a random type, then applying
// |func| on that type.  Used below with ValueGenerator to write values into a
// StructureSink, and with ValueChecker to validate values read from a
// StructureSource.
template <typename FuncT>
void IterateFuncWithType(FuncT func, int iterations) {
  PRNG type_generator;
  for (int i = 0; i < iterations; ++i) {
    // Pick a random type.
    static const int kNumTypes = 10;
    switch (std::uniform_int_distribution<int>(0, kNumTypes)(type_generator)) {
      case 0:
        func(TypeParam<int8>());
        break;
      case 1:
        func(TypeParam<uint8>());
        break;
      case 2:
        func(TypeParam<int16>());
        break;
      case 3:
        func(TypeParam<uint16>());
        break;
      case 4:
        func(TypeParam<int32>());
        break;
      case 5:
        func(TypeParam<uint32>());
        break;
      case 6:
        func(TypeParam<int64>());
        break;
      case 7:
        func(TypeParam<uint64>());
        break;
      case 8:
        func(TypeParam<float>());
        break;
      case 9:
        func(TypeParam<double>());
        break;
    }
  }
}

// Functor class to write random values into a StructureSink.
class ValueGenerator {
 public:
  ValueGenerator(PRNG* generator, StructureSink* sink)
      : generator_(generator), sink_(sink) {}

  // Writes a random value of the type encoded in the unnamed proxy parameter to
  // the sink.
  //
  // The proxy parameter exists only to provide the type to the operator, not
  // the value.
  template <typename T>
  void operator()(TypeParam<T>) {
    sink_->WritePod(Random<T>(generator_));
  }

 private:
  PRNG* const generator_;
  StructureSink* const sink_;
};

// Functor class to read values from a StructureSource, and compare them against
// a stream of values generated from a PRNG.
class ValueChecker {
 public:
  ValueChecker(PRNG* generator, StructureSource* source)
      : generator_(generator), source_(source), value_index_(0) {}

  // Reads a value of the type encoded in the unnamed proxy parameter from the
  // source, and compares it against a value from a known sequence of (random)
  // numbers.
  //
  // The proxy parameter exists only to provide the type to the operator, not
  // the value.
  template <typename T>
  void operator()(TypeParam<T>) {
    EXPECT_EQ(Random<T>(generator_), source_->ReadPod<T>())
        << "value_index_=" << value_index_
        << " of kFuzzTestValueCount=" << kFuzzTestValueCount;
    ++value_index_;
  }

 private:
  PRNG* const generator_;
  StructureSource* const source_;
  int value_index_;
};

// Fuzz tests a large stream of values into a StructureSink, and checks that the
// equivalent reading from a StructureSource will restore its original values.
TEST(StructuredIoTest, FuzzTestRandomValues) {
  std::string buf;
  base::StringByteSink byte_sink(&buf);
  StructureSink sink(&byte_sink);

  // I rolled a d20.  Guaranteed to be random.
  const int kSeed = 17;

  // Initialize |sink_generator| and |source_generator| with the same seed.  The
  // pseudo-random stream generated from the common initial state should be
  // identical, allowing us to verify the random values written into the sink
  // (generated using |sink_generator|) against the values read from the source
  // (compared against |source_generator|).
  PRNG sink_generator(kSeed);
  PRNG source_generator(kSeed);

  // Generate |kFuzzTestValueCount| values and write into |sink|.
  IterateFuncWithType(ValueGenerator(&sink_generator, &sink),
                      kFuzzTestValueCount);

  // Add some ad-hoc testing for non-POD types.
  for (int i = 0; i < 16; ++i) {
    sink.WriteVector(Random<ion::math::Vector2i>(&sink_generator));
    sink.WriteVector(Random<ion::math::Vector3f>(&sink_generator));
    sink.WritePoint(Random<ion::math::Point2i>(&sink_generator));
    sink.WritePoint(Random<ion::math::Point3f>(&sink_generator));
    sink.WriteColor(Random<Color4ui>(&sink_generator));
    sink.WriteMatrix(Random<ion::math::Matrix4d>(&sink_generator));
    sink.WriteString(Random<std::string>(&sink_generator));
    const std::string buffer = Random<std::string>(&sink_generator);
    sink.WriteBytes(buffer.data(), buffer.size());
  }

  // Set up a StructureSource to read from written data.
  ArrayByteSource byte_source(absl::string_view(buf.data(), buf.size()));
  StructureSource source(&byte_source);

  // Validate |kFuzzTestValueCount| values from |source|.
  IterateFuncWithType(ValueChecker(&source_generator, &source),
                      kFuzzTestValueCount);

  // Validate the ad-hoc non-POD data.
  for (int i = 0; i < 16; ++i) {
    EXPECT_EQ(Random<ion::math::Vector2i>(&source_generator),
              source.ReadVector<ion::math::Vector2i>());
    EXPECT_EQ(Random<ion::math::Vector3f>(&source_generator),
              source.ReadVector<ion::math::Vector3f>());
    EXPECT_EQ(Random<ion::math::Point2i>(&source_generator),
              source.ReadPoint<ion::math::Point2i>());
    EXPECT_EQ(Random<ion::math::Point3f>(&source_generator),
              source.ReadPoint<ion::math::Point3f>());
    EXPECT_EQ(Random<Color4ui>(&source_generator),
              source.ReadColor<Color4ui>());
    EXPECT_EQ(Random<ion::math::Matrix4d>(&source_generator),
              source.ReadMatrix<ion::math::Matrix4d>());
    EXPECT_EQ(Random<std::string>(&source_generator), source.ReadString());
    const std::string buffer = Random<std::string>(&source_generator);
    std::vector<char> read_buffer(buffer.size());
    source.ReadBytes(read_buffer.data(), read_buffer.size());
    EXPECT_EQ(0, memcmp(buffer.data(), read_buffer.data(), read_buffer.size()));
  }
}

}  // namespace
}  // namespace base
}  // namespace seurat
