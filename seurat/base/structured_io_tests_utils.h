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

#ifndef VR_SEURAT_BASE_STRUCTURED_IO_TESTS_UTILS_H_
#define VR_SEURAT_BASE_STRUCTURED_IO_TESTS_UTILS_H_

#include <array>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <type_traits>

#include "ion/math/matrix.h"
#include "ion/math/vector.h"
#include "seurat/base/structured_io.h"
#include "seurat/base/util.h"

// This file contains utilities for Structured IO tests in seurat/base.

namespace seurat {
namespace base {

// These templated overloads of Random<> generate random numbers for testing,
// with a bias towards "special" values (min, max, zero, infinity, etc.).
// Overloads are also provided which generate ion::math::Vector,
// ion::math::Matrix, etc. types.  The functions take as an argument a random
// number generator as the source of randomness.

// Generate an integral type random.
template <typename T, typename RNG>
typename std::enable_if<std::is_integral<T>::value, T>::type Random(
    RNG* generator) {
  static const std::array<T, 5> special_values = {{
      std::numeric_limits<T>::min(),
      static_cast<T>(-1),
      static_cast<T>(0),
      static_cast<T>(1),
      std::numeric_limits<T>::max(),
  }};
  if (std::uniform_int_distribution<int>(0, 3)(*generator) < 1) {
    // Return a random special value some fraction of the time.
    return special_values[std::uniform_int_distribution<size_t>(
        0, special_values.size() - 1)(*generator)];
  } else {
    // Return just a random value.
    //
    // Note this uses the DistributionTrait type to replace T with an
    // appropriate type for the operating precision of the distribution.
    using DistributionType =
        typename base::DistributionTrait<T>::DistributionOperatingType;
    return static_cast<T>(std::uniform_int_distribution<DistributionType>(
        std::numeric_limits<T>::min(),
        std::numeric_limits<T>::max())(*generator));
  }
}

// Generate a floating-point type random.
template <typename T, typename RNG>
typename std::enable_if<std::is_floating_point<T>::value, T>::type Random(
    RNG* generator) {
  static const std::array<T, 10> special_values = {{
      std::numeric_limits<T>::min(),
      static_cast<T>(-1.0),
      static_cast<T>(-0.0),
      static_cast<T>(0.0),
      static_cast<T>(1.0),
      std::numeric_limits<T>::max(),
      std::numeric_limits<T>::epsilon(),
      std::numeric_limits<T>::infinity(),
      -std::numeric_limits<T>::infinity(),
      std::numeric_limits<T>::denorm_min(),
  }};
  if (std::uniform_int_distribution<int>(0, 3)(*generator) < 1) {
    // Return a random special value some fraction of the time.
    return special_values[std::uniform_int_distribution<size_t>(
        0, special_values.size() - 1)(*generator)];
  } else {
    // Return just a random value.
    return std::uniform_real_distribution<T>(
        std::numeric_limits<T>::min(),
        std::numeric_limits<T>::max())(*generator);
  }
}

// Generate an std::string random.
template <typename T, typename RNG>
typename std::enable_if<
    std::is_same<typename std::decay<T>::type, std::string>::value, T>::type
Random(RNG* generator) {
  size_t length = std::uniform_int_distribution<size_t>(0, 1024)((*generator));
  std::unique_ptr<char[]> buffer(new char[length]);
  for (size_t i = 0; i < length; ++i) {
    buffer[i] = Random<char>(generator);
  }
  return std::string(buffer.get(), length);
}

// Generate an ion::math::VectorBase<> type random (Point, Vector).
template <typename T, typename RNG>
typename std::enable_if<
    std::is_convertible<typename std::decay<T>::type,
                        typename ion::math::VectorBase<
                            T::kDimension, typename T::ValueType>>::value,
    T>::type
Random(RNG* generator) {
  T ret;
  for (int i = 0; i < T::kDimension; ++i) {
    ret[i] = Random<typename T::ValueType>(generator);
  }
  return ret;
}

// Generate an ion::math::Matrix type random.
template <typename T, typename RNG>
typename std::enable_if<
    std::is_same<typename std::decay<T>::type,
                 typename ion::math::Matrix<T::kDimension,
                                            typename T::ValueType>>::value,
    T>::type
Random(RNG* generator) {
  T ret;
  for (int i = 0; i < T::kDimension; ++i) {
    for (int j = 0; j < T::kDimension; ++j) {
      ret[i][j] = Random<typename T::ValueType>(generator);
    }
  }
  return ret;
}

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_STRUCTURED_IO_TESTS_UTILS_H_
