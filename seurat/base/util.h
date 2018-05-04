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

#ifndef VR_SEURAT_BASE_UTIL_H_
#define VR_SEURAT_BASE_UTIL_H_

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>

namespace seurat {
namespace base {

namespace util_internal {

template <typename ToType, typename FromType, int Index>
struct Converter {
  template <typename... Args>
  static ToType Convert(const FromType& f, Args&&... args) {
    return Converter<ToType, FromType, Index - 1>::Convert(
        f, f[Index - 1], std::forward<Args>(args)...);
  }
};

template <typename ToType, typename FromType>
struct Converter<ToType, FromType, 0> {
  template <typename... Args>
  static ToType Convert(const FromType& f, Args&&... args) {
    return ToType(std::forward<Args>(args)...);
  }
};

}  // namespace util_internal

// Converts from one vector type to another, of the same dimension.  The source
// vector type must have a |kDimension| static member.
template <typename ToType, typename FromType>
ToType ConvertVector(const FromType& from) {
  return util_internal::Converter<ToType, FromType,
                                  FromType::kDimension>::Convert(from);
}

// DistributionTrait<T> defines types used for operating space for random number
// distributions that output other, typically lower-range, types.
template <typename T>
struct DistributionTrait {
  // By default, distributions should use the same type as the output value.
  // However, for byte-size types, we override the default choice with int.
  // This is because uniform_int_distribution<uint8> is non-standard, according
  // to [rand.req.genl] (26.5.1.1) and MSVC obeys this rule but GCC does not:
  //
  // "Throughout this subclause 26.5, the effect of instantiating a template:
  // ... that has a template type parameter named IntType is undefined unless
  // the corresponding template argument is cv-unqualified and is one of short,
  // int, long, long long, unsigned short, unsigned int, unsigned long, or
  // unsigned long long."
  using DistributionOperatingType =
      typename std::conditional<sizeof(T) >= sizeof(std::int16_t), T,
                                std::int16_t>::type;
};

// Retrieves the operating type given an integral type.
template <typename T>
using DistributionOpType =
    typename DistributionTrait<T>::DistributionOperatingType;

// Rounds the input |number| up to a multiple of |n|.
template <typename T>
inline T RoundModN(T number, T n) {
  number += (n - 1);
  number -= (number % n);
  return number;
}

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_UTIL_H_
