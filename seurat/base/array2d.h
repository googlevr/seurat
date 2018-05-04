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

#ifndef VR_SEURAT_BASE_ARRAY2D_H_
#define VR_SEURAT_BASE_ARRAY2D_H_

#include <algorithm>
#include <type_traits>
#include <vector>

#include "ion/math/vector.h"

namespace seurat {
namespace base {

// A 2D array of type T, e.g. seurat::base::Color3f. Array2D
// data has an orientation, bottom-up, Cartesian coordinates. Increasing memory
// addresses therefore correspond to increasing Y (or V) coordinates.
template <typename T>
class Array2D {
 public:
  using ElementType = T;

  // Creates an empty 0x0 array.
  Array2D() : width_(0), height_(0) {}

  // Creates an array of the specified size with default-constructed values.
  Array2D(int width, int height)
      : data_(width * height), width_(width), height_(height) {}

  // Creates an array of the specified |size| with default-constructed values.
  explicit Array2D(const ion::math::Vector2i& size)
      : Array2D(size[0], size[1]) {}

  // Creates an array of the specified size initialized with the given |value|.
  Array2D(int width, int height, const ElementType& value)
      : data_(width * height, ElementData{value}),
        width_(width),
        height_(height) {}

  // Creates an array of the specified |size| initialized with the given
  // |value|.
  Array2D(const ion::math::Vector2i& size, const ElementType& value)
      : Array2D(size[0], size[1], value) {}

  // Creates an array by copying the contents of the given lvalue |rhs|.
  Array2D(const Array2D<T>& rhs) : width_(0), height_(0) { operator=(rhs); }

  // Creates an array by moving the contents of the given rvalue |rhs|, leaving
  // |rhs| empty.
  Array2D(Array2D<T>&& rhs) : width_(0), height_(0) {
    operator=(std::move(rhs));
  }

  // Resizes the array. The elements are initialized with the default
  // constructor.
  void Resize(int width, int height) {
    data_.resize(width * height);
    width_ = width;
    height_ = height;
  }

  // Resizes the array. The elements are initialized with the default
  // constructor.
  void Resize(const ion::math::Vector2i& size) { Resize(size[0], size[1]); }

  // Resizes the array and initialize all elements with |value|.
  void Resize(int width, int height, const ElementType& value) {
    data_.resize(width * height, ElementData{value});
    width_ = width;
    height_ = height;
  }

  // Resizes the array and initialize all elements with |value|.
  void Resize(const ion::math::Vector2i& size, const ElementType& value) {
    Resize(size[0], size[1], value);
  }

  // Copies the contents of the given lvalue |rhs| into this.
  Array2D<T>& operator=(const Array2D<T>& rhs) {
    data_ = rhs.data_;
    width_ = rhs.width_;
    height_ = rhs.height_;
    return *this;
  }

  // Moves the contents of the given rvalue |rhs| into this, leaving |rhs|
  // empty.
  Array2D<T>& operator=(Array2D<T>&& rhs) {
    data_ = std::move(rhs.data_);
    width_ = rhs.width_;
    height_ = rhs.height_;
    rhs.width_ = 0;
    rhs.height_ = 0;
    return *this;
  }

  // Compares for equality, in size and element values.
  bool operator==(const Array2D<T>& other) const;

  // Compares for inequality, in size and element values.
  bool operator!=(const Array2D<T>& other) const;

  // Sets all elements to the specified value.
  void Fill(const ElementType& value);

  // Returns the width of the array.
  int Width() const { return width_; }

  // Returns the height of the array.
  int Height() const { return height_; }

  // Returns the size of the array.
  ion::math::Vector2i GetSize() const { return {width_, height_}; }

  // Returns true iff the |p| is inside the Array2D.
  bool IsInside(const ion::math::Point2i& p) const {
    return 0 <= p[0] && p[0] < width_ && 0 <= p[1] && p[1] < height_;
  }

  // Returns a const reference to the element at (x, y).
  const ElementType& At(int x, int y) const;

  // Returns a mutable reference to the element at (x, y).
  ElementType& At(int x, int y);

  // Returns a const reference to the element at coordinates (p[0], p[1]).
  const ElementType& At(const ion::math::Point2i& p) const {
    return At(p[0], p[1]);
  }

  // Returns a mutable reference to the element at coordinates (p[0], p[1]).
  ElementType& At(const ion::math::Point2i& p) { return At(p[0], p[1]); }

  // Returns a const pointer to raw element data.
  //
  // Only enabled if T != bool.
  const ElementType* Data() const { return data_.data(); }

  // Returns a pointer to raw element data.
  //
  // Only enabled if T != bool.
  ElementType* Data() { return data_.data(); }

  // The total number of pixels.
  size_t size() const { return data_.size(); }

  // STL-style begin(), end(), and data() to support Span & foreach loops.
  //
  // Only enabled if T != bool.
  const ElementType* data() const { return data_.data(); }
  ElementType* data() { return data_.data(); }
  typename std::vector<ElementType>::const_iterator begin() const {
    return data_.begin();
  }
  typename std::vector<ElementType>::const_iterator end() const {
    return data_.end();
  }
  typename std::vector<ElementType>::iterator begin() { return data_.begin(); }
  typename std::vector<ElementType>::iterator end() { return data_.end(); }

 private:
  // To support ElementType==bool, wrap ElementType with a struct to force
  // std::vector to not specialize on bool.
  struct ElementWrapper {
    ElementType value;

    bool operator==(const ElementWrapper& other) const {
      return value == other.value;
    }
  };

  // Only use the wrapper to workaround ElementType==bool.
  using ElementData =
      typename std::conditional<std::is_same<ElementType, bool>::value,
                                ElementWrapper, ElementType>::type;

  // Convert ElementData into ElementType via method overloading.
  //
  // In the case of the wrapper, we must unwrap it.
  static ElementType* Unwrap(ElementWrapper* wrapper) {
    return &wrapper->value;
  }
  static const ElementType* Unwrap(const ElementWrapper* wrapper) {
    return &wrapper->value;
  }
  // If there is no wrapper, then do nothing.
  static ElementType* Unwrap(ElementType* unwrapped) { return unwrapped; }
  static const ElementType* Unwrap(const ElementType* unwrapped) {
    return unwrapped;
  }

  // Data in scanline order.
  std::vector<ElementData> data_;

  // The width of the array.
  int width_;

  // The height of the array.
  int height_;
};

template <typename T>
bool Array2D<T>::operator==(const Array2D<T>& other) const {
  if (width_ != other.width_) {
    return false;
  }
  if (height_ != other.height_) {
    return false;
  }
  return data_ == other.data_;
}

template <typename T>
bool Array2D<T>::operator!=(const Array2D<T>& other) const {
  return !operator==(other);
}

template <typename T>
void Array2D<T>::Fill(const ElementType& value) {
  std::fill(data_.begin(), data_.end(), ElementData{value});
}

template <typename T>
const typename Array2D<T>::ElementType& Array2D<T>::At(int x, int y) const {
  DCHECK(IsInside({x, y})) << ion::math::Point2i(x, y);
  return *Unwrap(&data_[width_ * y + x]);
}

template <typename T>
typename Array2D<T>::ElementType& Array2D<T>::At(int x, int y) {
  DCHECK(IsInside({x, y})) << ion::math::Point2i(x, y);
  return *Unwrap(&data_[width_ * y + x]);
}

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_ARRAY2D_H_
