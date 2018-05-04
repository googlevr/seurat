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

#ifndef VR_SEURAT_BASE_ARRAY2D_VIEW_H_
#define VR_SEURAT_BASE_ARRAY2D_VIEW_H_

#include <algorithm>

#include "ion/math/vector.h"
#include "seurat/base/array2d.h"

namespace seurat {
namespace base {

// A mutable view into a dense 2D array.
template <typename T>
class MutableArray2DView {
 public:
  using ElementType = T;

  // Construct an empty MutableArray2DView.
  MutableArray2DView() : data_(nullptr), width_(0), height_(0), stride_(0) {}

  // Construct an MutableArray2DView which indexes into the |data| according
  // to the specified row |stride|.
  MutableArray2DView(T* data, int width, int height, int stride)
      : data_(data), width_(width), height_(height), stride_(stride) {}

  // Construct an MutableArray2DView which indexes into the |data| assuming no
  // padding.
  MutableArray2DView(T* data, const ion::math::Vector2i& size)
      : MutableArray2DView(data, size[0], size[1], size[0]) {}

  // Crop an existing MutableArray2DView to the subregion starting at index
  // (x, y) with the specified dimensions.
  MutableArray2DView(MutableArray2DView<T> original, int x, int y, int width,
                     int height);

  // Construct an MutableArray2DView of a subregion of an Array2D.
  MutableArray2DView(Array2D<T>* array2d, int x, int y, int width, int height)
      : MutableArray2DView(MutableArray2DView(array2d), x, y, width, height) {}

  // Construct an MutableArray2DView of an Array2D.
  explicit MutableArray2DView(Array2D<T>* array2d)
      : MutableArray2DView(array2d->Data(), array2d->Width(), array2d->Height(),
                           array2d->Width()) {}

  // Compares for equality, in size and element values.
  bool operator==(const MutableArray2DView<T>& other) const;

  // Compares for inequality, in size and element values.
  bool operator!=(const MutableArray2DView<T>& other) const;

  // Sets all elements to the specified value.
  void Fill(const ElementType& value);

  // Returns the width of the array.
  int Width() const { return width_; }

  // Returns the height of the array.
  int Height() const { return height_; }

  // Returns the size of the array.
  ion::math::Vector2i GetSize() const { return {width_, height_}; }

  // Returns true iff the |p| is inside the MutableArray2DView.
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

  // Returns the stride of data, e.g. for converting view types.
  int Stride() const { return stride_; }

 private:
  // A pointer to the element at index (0, 0).
  T* data_;
  // The width of the logical 2D array represented by this view.
  int width_;
  // The height of the logical 2D array represented by this view.
  int height_;
  // The number of indices in the data_ between rows in the view.
  //
  // In other words, if the value of index (x, y) has address I, then I +
  // stride_ is the address of (x, y + 1).
  int stride_;
};

// An immutable view into a dense 2D array.
template <typename T>
class Array2DView : public MutableArray2DView<const T> {
 public:
  // Construct an empty Array2DView.
  Array2DView() : MutableArray2DView<const T>() {}

  // Construct an Array2DView which wraps the image |data| with the specified
  // dimensiona and stride..
  Array2DView(T const* data, int width, int height, int stride)
      : MutableArray2DView<const T>(data, width, height, stride) {}
  // Construct an Array2DView which wraps the image |data| with the specified
  // dimensions and stride.
  Array2DView(T const* data, const ion::math::Vector2i& size)
      : Array2DView(data, size[0], size[1], size[0]) {}
  // Construct a view into a subregion of the |original| view.
  Array2DView(Array2DView<T> original, int x, int y, int width, int height)
      : MutableArray2DView<const T>(original, x, y, width, height) {}
  // Construct a view into a subregion of the |original| array2d.
  Array2DView(const Array2D<T>* array2d, int x, int y, int width, int height)
      : MutableArray2DView<const T>(Array2DView(array2d), x, y, width, height) {
  }
  // Construct a view wrapping the original |array2d|.
  explicit Array2DView(const Array2D<T>* array2d)
      : MutableArray2DView<const T>(array2d->Data(), array2d->Width(),
                                    array2d->Height(), array2d->Width()) {}
  // Converts mutable view |original| to an identical immutable view.
  explicit Array2DView(const MutableArray2DView<T>& original)
      : MutableArray2DView<const T>(&original.At(0, 0), original.Width(),
                                    original.Height(), original.Stride()) {}
};

template <typename T>
MutableArray2DView<T>::MutableArray2DView(MutableArray2DView<T> original, int x,
                                          int y, int width, int height) {
  CHECK(original.IsInside({x, y})) << "Invalid (x, y) coordinates for crop";
  CHECK_GE(width, 0) << "Invalid width";
  CHECK_GE(height, 0) << "Invalid height";
  CHECK_LE(x + width, original.Width()) << "Invalid crop width";
  CHECK_LE(y + height, original.Height()) << "Invalid crop height";

  data_ = original.data_ + y * original.stride_ + x;
  width_ = width;
  height_ = height;
  stride_ = original.stride_;
}

template <typename T>
bool MutableArray2DView<T>::operator==(
    const MutableArray2DView<T>& other) const {
  if (width_ != other.width_) {
    return false;
  }
  if (height_ != other.height_) {
    return false;
  }
  for (int y = 0; y < height_; ++y) {
    if (!std::equal(&data_[stride_ * y], &data_[stride_ * y + width_],
                    &other.data_[other.stride_ * y])) {
      return false;
    }
  }
  return true;
}

template <typename T>
bool MutableArray2DView<T>::operator!=(
    const MutableArray2DView<T>& other) const {
  return !operator==(other);
}

template <typename T>
void MutableArray2DView<T>::Fill(const ElementType& value) {
  for (int y = 0; y < height_; ++y) {
    std::fill(&data_[stride_ * y], &data_[stride_ * y + width_], value);
  }
}

template <typename T>
const typename MutableArray2DView<T>::ElementType& MutableArray2DView<T>::At(
    int x, int y) const {
  DCHECK(IsInside({x, y})) << ion::math::Point2i(x, y);
  return data_[stride_ * y + x];
}

template <typename T>
typename MutableArray2DView<T>::ElementType& MutableArray2DView<T>::At(int x,
                                                                       int y) {
  DCHECK(IsInside({x, y})) << ion::math::Point2i(x, y);
  return data_[stride_ * y + x];
}

}  // namespace base
}  // namespace seurat

#endif  // VR_SEURAT_BASE_ARRAY2D_VIEW_H_
