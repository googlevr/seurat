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

#ifndef VR_SEURAT_BASE_STRUCTURED_IO_H_
#define VR_SEURAT_BASE_STRUCTURED_IO_H_

#include <cstdint>
#include <numeric>
#include <string>
#include <type_traits>
#include <vector>

#include "ion/gfx/image.h"
#include "ion/math/matrix.h"
#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "seurat/base/bytestream.h"
#include "seurat/base/color.h"

namespace seurat {
namespace base {

// Supports WriteObjectsToSink & ReadObjectsFromSink prototypes.
class StructureSink;
class StructureSource;

namespace structured_io_internal {

// WriteObjectsToSink assists in serializing user-defined types. The internal
// definition allows any plain-old-data type. StructureSink's WritePod methods
// support any type that wants to "opt in" to binary serialization by providing
// a WriteObjectsToSink implementation. Additionally, the user defined type may
// provide a WriteObjectsToSink implementation with member-by-member
// serialization as needed for polymorphic types. StructureSink finds such
// implementations via argument-dependent lookup, analogous to std::swap's
// implementation. This version handles global types like int and char, in
// addition to aggregate types.
//
// N.B. if lookup fails for an overload of one of these helpers, ensure the
// function is implemented in the same namespace as the type. For example,
// adding an implementation for a naked type from the std namespace will fail,
// because it would need to be inside std, but the C++ standard rules ban most
// such declarations as undefined behavior. One workaround is to use a wrapper
// type instead.
template <typename T>
inline typename std::enable_if<std::is_pod<T>::value, void>::type
WriteObjectsToSink(const T* objects, size_t count, StructureSink* source);

// Assists in deserializing user-defined types, like WriteObjectsToSink. See
// WriteObjectsToSink for important details.
template <typename T>
inline typename std::enable_if<std::is_pod<T>::value, void>::type
ReadObjectsFromSource(StructureSource* source, T* objects, size_t count);

}  // namespace structured_io_internal

// Stores structured data - in contrast to raw byte streams - e.g. for
// serializing objects.
//
// Data must be read via StructureSink in the order in which it was
// written.
class StructureSink {
 public:
  // Constructs a sink which writes to the given ByteSink.
  explicit StructureSink(ByteSink* sink) : sink_(sink) {}

  void WriteBytes(const char* bytes, size_t count);

  // Low-level IO allows writing an array of objects of any type to the sink.
  // This is only intended to provide type safety to WriteObjectsToSink
  // implementations. It's public to allow WriteObjectsToSink implementations
  // for types not known to StructureSink to use it.
  template <typename T>
  void WriteObjects(const T* objects, int64_t count) {
    WriteBytes(reinterpret_cast<const char*>(objects), sizeof(T) * count);
  }

  void WriteString(const std::string& str);

  template <typename T>
  typename std::enable_if<std::is_pod<T>::value, void>::type WritePod(
      const T& value) {
    static_assert(!std::is_pointer<T>::value,
                  "Correct WritePod takes const& params, not pointers. There "
                  "may be an accidental address-of (&) operator on the "
                  "parameter.");
    WriteBytes(reinterpret_cast<const char*>(&value), sizeof(T));
  }

  // Write a fixed-size array of |Dimension| pod values from |elements|,
  // where the array's type implies |Dimension| at compile time.
  // This overload only handles 1D arrays, and the enable_if enforces that.
  template <int Dimension, typename T>
  typename std::enable_if<!std::is_array<T>::value, void>::type WritePodArray(
      const T (&elements)[Dimension]) {
    static_assert(std::is_pod<T>::value && Dimension > 0,
                  "WritePodArray requires a plain-old-data element type");
    WritePodArray(elements, Dimension);
  }

  // Write a multidimensional array |elements| of pod values, where the array's
  // type implies the outer and inner sizes at compile time.
  //
  // Note: The inner array may in turn be multidimensional, and this overload
  // recursively removes outer arrays by writing the set of inner arrays with
  // WritePodArray.
  template <int DimensionO, int DimensionI, typename T>
  void WritePodArray(const T (&elements)[DimensionO][DimensionI]) {
    static_assert(std::is_pod<T>::value,
                  "WritePodArray requires a plain-old-data element type");
    // Take off the outer array.
    for (size_t outer_index = 0; outer_index < DimensionO; ++outer_index)
      WritePodArray(elements[outer_index]);
  }

  // Write a variable-size array of |dimension| pod data items from |elements|,
  // where |dimension| is determined at run time by the parameter.
  template <typename T>
  void WritePodArray(const T* elements, size_t dimension) {
    static_assert(std::is_pod<T>::value,
                  "WritePodArray requires a plain-old-data element type");
    for (size_t element_index = 0; element_index < dimension; ++element_index)
      WritePod<T>(elements[element_index]);
  }

  // Write a variable-sized vector of pod values into |elements|.
  template <typename T>
  void WritePodArray(const std::vector<T>& elements) {
    CHECK_LE(elements.size(),
             static_cast<size_t>(std::numeric_limits<int64_t>::max()));
    WritePod<int64_t>(static_cast<int64_t>(elements.size()));
    // Import this name to use default implementation, but allow argument-
    // dependent lookup to find external implementations in other namespaces,
    // like std::swap.
    //
    // The default implementation supports basic global types like int, and it
    // must be made visible because it cannot be found through ADL for global
    // types as they are unnamespaced.
    using structured_io_internal::WriteObjectsToSink;
    WriteObjectsToSink(elements.data(), elements.size(), this);
  }

  template <int Dimension, typename T>
  void WriteVector(const ion::math::Vector<Dimension, T>& vec) {
    WritePodArray(vec.Data(), Dimension);
  }

  template <int Dimension, typename T>
  void WritePoint(const ion::math::Point<Dimension, T>& point) {
    WritePodArray(point.Data(), Dimension);
  }

  template <int Dimension, typename T>
  void WriteColor(const Color<Dimension, T>& color) {
    WritePodArray(color.Data(), Dimension);
  }

  template <int Dimension, typename T>
  void WriteMatrix(const ion::math::Matrix<Dimension, T>& mat) {
    WritePodArray(mat.Data(), Dimension * Dimension);
  }

  void WriteImage(const ion::gfx::ImagePtr image);

 private:
  ByteSink* sink_;

  DISALLOW_COPY_AND_ASSIGN(StructureSink);
};

// Used for deserializing structured data, e.g. objects, written via
// StructureSink.
class StructureSource {
 public:
  // Constructs a sink which reads from the given ByteSink.
  explicit StructureSource(ByteSource* source) : source_(source) {}

  void ReadBytes(char* bytes_out, size_t count);

  // Low-level IO allows reading an array of objects of any type from the
  // source. This is only intended to provide type safety to ReadObjectsFromSink
  // implementations. It's public to allow ReadObjectsFromSink implementations
  // for types not known to StructureSource to use it.
  template <typename T>
  void ReadObjects(T* objects, int64_t count) {
    ReadBytes(reinterpret_cast<char*>(objects), sizeof(T) * count);
  }

  std::string ReadString();

  // Template magic.  Call with:
  // POD_TYPE pod = source.ReadPod<POD_TYPE>();
  template <typename T>
  typename std::enable_if<std::is_pod<typename std::decay<T>::type>::value,
                          typename std::decay<T>::type>::type
  ReadPod() {
    typename std::decay<T>::type value;
    ReadBytes(reinterpret_cast<char*>(&value), sizeof(T));
    return value;
  }

  // Read a fixed-size array of |Dimension| pod values from |elements|,
  // where the array's type implies |Dimension| at compile time.
  // This overload only handles 1D arrays, and the enable_if enforces that.
  template <int Dimension, typename T>
  typename std::enable_if<!std::is_array<T>::value, void>::type ReadPodArray(
      T (&elements)[Dimension]) {
    static_assert(std::is_pod<T>::value && Dimension > 0,
                  "ReadPodArray requires a plain-old-data element type");
    ReadPodArray(elements, Dimension);
  }

  // Read a multidimensional array |elements| of pod values, where the array's
  // type implies the outer and inner sizes at compile time.
  //
  // Note: The inner array may in turn be multidimensional, and this overload
  // recursively removes outer arrays by reading the set of inner arrays with
  // ReadPodArray.
  template <int DimensionO, int DimensionI, typename T>
  void ReadPodArray(T (&elements)[DimensionO][DimensionI]) {
    static_assert(std::is_pod<T>::value,
                  "ReadPodArray requires a plain-old-data element type");
    // Take off the outer array.
    for (size_t outer_index = 0; outer_index < DimensionO; ++outer_index)
      ReadPodArray(elements[outer_index]);
  }

  // Read a variable-sized array of |dimension| pod values into |elements|,
  // where the |dimension| parameter specifies the count of elements.
  template <typename T>
  void ReadPodArray(T* elements, size_t dimension) {
    static_assert(std::is_pod<T>::value,
                  "ReadPodArray requires a plain-old-data element type");
    for (size_t element_index = 0; element_index < dimension; ++element_index)
      elements[element_index] = ReadPod<T>();
  }

  // Read a variable-sized vector of pod values into |elements|.
  template <typename T>
  void ReadPodArray(std::vector<T>* elements) {
    const int64_t element_count = ReadPod<int64_t>();
    CHECK_GE(element_count, 0);
    elements->resize(element_count);
    // Import this name to use default implementation, but allow argument
    // dependent lookup to find external implementations in other namespaces,
    // like std::swap.
    //
    // The default implementation supports basic global types like int, and it
    // must be made visible because it cannot be found through ADL for global
    // types as they are unnamespaced.
    using structured_io_internal::ReadObjectsFromSource;
    ReadObjectsFromSource(this, elements->data(), element_count);
  }

  // Template magic.  Call with:
  // VectorXX vector = source.ReadVector<VectorXX>();
  template <typename T>
  typename std::enable_if<
      std::is_same<typename std::decay<T>::type,
                   typename ion::math::Vector<T::kDimension,
                                              typename T::ValueType>>::value,
      typename std::decay<T>::type>::type
  ReadVector() {
    typename std::decay<T>::type vec;
    ReadPodArray(vec.Data(), T::kDimension);
    return vec;
  }

  // Template magic.  Call with:
  // PointXX point = source.ReadPoint<PointXX>();
  template <typename T>
  typename std::enable_if<
      std::is_same<typename std::decay<T>::type,
                   typename ion::math::Point<T::kDimension,
                                             typename T::ValueType>>::value,
      typename std::decay<T>::type>::type
  ReadPoint() {
    typename std::decay<T>::type point;
    ReadPodArray(point.Data(), T::kDimension);
    return point;
  }

  // Template magic.  Call with:
  // PointXX point = source.ReadPoint<PointXX>();
  template <typename T>
  typename std::enable_if<
      std::is_same<
          typename std::decay<T>::type,
          typename base::Color<T::kDimension, typename T::ValueType>>::value,
      typename std::decay<T>::type>::type
  ReadColor() {
    typename std::decay<T>::type color;
    ReadPodArray(color.Data(), T::kDimension);
    return color;
  }

  // Template magic.  Call with:
  // MatrixXX matrix = source.ReadMatrix<MatrixXX>();
  template <typename T>
  typename std::enable_if<
      std::is_convertible<typename std::decay<T>::type,
                          typename ion::math::Matrix<
                              T::kDimension, typename T::ValueType>>::value,
      typename std::decay<T>::type>::type
  ReadMatrix() {
    typename std::decay<T>::type mat;
    ReadPodArray(mat.Data(), T::kDimension * T::kDimension);
    return mat;
  }

  ion::gfx::ImagePtr ReadImage();

 private:
  ByteSource* source_;

  DISALLOW_COPY_AND_ASSIGN(StructureSource);
};

namespace structured_io_internal {

template <typename T>
typename std::enable_if<std::is_pod<T>::value, void>::type WriteObjectsToSink(
    const T* objects, size_t count, StructureSink* sink) {
  sink->WriteObjects(objects, count);
}

template <typename T>
typename std::enable_if<std::is_pod<T>::value, void>::type
ReadObjectsFromSource(StructureSource* source, T* objects, size_t count) {
  source->ReadObjects(objects, count);
}

}  // namespace structured_io_internal

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_STRUCTURED_IO_H_
