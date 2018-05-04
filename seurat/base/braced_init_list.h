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

#ifndef VR_SEURAT_BASE_BRACED_INIT_LIST_H_
#define VR_SEURAT_BASE_BRACED_INIT_LIST_H_

#include <array>
#include <sstream>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <vector>

#include "ion/math/vector.h"
#include "seurat/base/array2d.h"
#include "seurat/base/color.h"

namespace seurat {
namespace base {

// ToBracedInitList() converts Points, Vectors, Colors, Polygons as well as
// possibly nested std::arrays or std::vectors to strings in the
// brace-init-list// format, for exporting data to Mathematica or C/C++
// programs.

// Specialization for integer and floating point types.
template <typename T, typename std::enable_if<
                          std::is_arithmetic<T>::value>::type* = nullptr>
std::string ToBracedInitList(T x) {
  std::stringstream buf;
  buf << +x;
  return buf.str();
}

// Specialization for classes derived from ion::math::VectorBase<D,T>, e.g.
// Vector2i, Point3f, Color4d.
template <int D, typename T>
std::string ToBracedInitList(const ion::math::VectorBase<D, T>& x) {
  std::stringstream buf;
  buf << "{";
  for (int i = 0; i < D; ++i) {
    if (0 < i) buf << " ";
    buf << +x[i];
    if (i < D - 1) buf << ",";
  }
  buf << "}";
  return buf.str();
}

// Specialization for Array2D and Image types.
template <typename T>
std::string ToBracedInitList(Array2D<T> array) {
  std::ostringstream buf;
  buf << "{";
  for (int y = 0; y < array.Height(); ++y) {
    buf << "{";
    for (int x = 0; x < array.Width(); ++x) {
      if (0 < x) buf << " ";
      buf << ToBracedInitList(array.At(x, y));
      if (x < array.Width() - 1) buf << ",";
    }
    buf << "}";
    if (y < array.Height() - 1) buf << ",\n";
  }
  buf << "}";
  return buf.str();
}

// Specialization for std-container-like classes.
template <typename T, typename T::iterator* = nullptr>
std::string ToBracedInitList(T x) {
  std::ostringstream buf;
  buf << "{";
  typename T::iterator last = x.end() - 1;
  for (typename T::iterator it = x.begin(); it != x.end(); ++it) {
    if (it != x.begin()) buf << " ";
    buf << ToBracedInitList(*it);
    if (it != last) buf << ",";
  }
  buf << "}";
  return buf.str();
}

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_BRACED_INIT_LIST_H_
