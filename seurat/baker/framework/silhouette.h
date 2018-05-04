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

#ifndef VR_SEURAT_BAKER_FRAMEWORK_SILHOUETTE_H_
#define VR_SEURAT_BAKER_FRAMEWORK_SILHOUETTE_H_

#include <memory>
#include <vector>

#include "ion/math/vector.h"
#include "seurat/geometry/kdtree.h"

namespace seurat {
namespace baker {

// A 2D Silhouette, implied by nearest-neighbors from a set of solid and
// freespace samples.
//
// In other words, a point is a 'solid' part of the silhouette iff it is closer
// to a solid sample than a freespace sample.
//
// Samples are typically in the range [0, 1]^2.
class ImplicitSilhouette {
 public:
  ImplicitSilhouette(std::vector<ion::math::Point2f> solid,
                     std::vector<ion::math::Point2f> freespace);

  // Returns whether the implicit silhouette is "solid" at the given |point|.
  bool IsSolidAtPoint(const ion::math::Point2f& point) const;

 private:
  // All solid samples.
  //
  // Must be const because solid_tree_ points into it.
  const std::vector<ion::math::Point2f> solid_;

  // All freespace samples.
  //
  // Must be const because freespace_tree_ points into it.
  const std::vector<ion::math::Point2f> freespace_;

  // An acceleration structure for querying solid_.
  const geometry::KdTree<2> solid_tree_;

  // An acceleration structure for querying freespace_.
  const geometry::KdTree<2> freespace_tree_;
};

// Builds an ImplicitSilhouette.
//
// For the sake of limiting memory consumption, different implementations
// have different tradeoffs (i.e. by discarding/merging samples).
class SilhouetteBuffer {
 public:
  virtual ~SilhouetteBuffer() = default;

  // Adds a solid sample to the buffer.
  virtual void AddSolidSample(const ion::math::Point2f& sample) = 0;

  // Adds a freespace sample to the buffer.
  virtual void AddFreespaceSample(const ion::math::Point2f& sample) = 0;

  // Returns the ImplicitSilhouette resulting from all previously added (and
  // possibly merged/discarded) samples, including those added before the last
  // Resolve().
  virtual std::unique_ptr<ImplicitSilhouette> Resolve() const = 0;
};

// The simplest SilhouetteBuffer.
//
// Retains all solid & freespace samples & adds them to the resulting
// ImplicitSilhouette.
class SimpleSilhouetteBuffer : public SilhouetteBuffer {
 public:
  SimpleSilhouetteBuffer() = default;
  ~SimpleSilhouetteBuffer() override = default;

  // SilhouetteBuffer implementation.
  void AddSolidSample(const ion::math::Point2f& sample) override {
    solid_.push_back(sample);
  }
  void AddFreespaceSample(const ion::math::Point2f& sample) override {
    freespace_.push_back(sample);
  }
  std::unique_ptr<ImplicitSilhouette> Resolve() const override {
    return std::unique_ptr<ImplicitSilhouette>(
        new ImplicitSilhouette(solid_, freespace_));
  }

 private:
  std::vector<ion::math::Point2f> solid_;
  std::vector<ion::math::Point2f> freespace_;
};

// A limited-memory SilhouetteBuffer.
//
// Retains at-most-one quantized solid & freespace sample per pixel over a
// predefined resolution grid.
//
// Note:  Only samples within the unit square [0.0, 1.0)^2 are retained.
class CompactSilhouetteBuffer : public SilhouetteBuffer {
 public:
  explicit CompactSilhouetteBuffer(const ion::math::Vector2i& resolution)
      : resolution_(resolution),
        solid_(resolution[0] * resolution[1]),
        freespace_(resolution[0] * resolution[1]) {}
  ~CompactSilhouetteBuffer() override = default;

  // SilhouetteBuffer implementation.
  void AddSolidSample(const ion::math::Point2f& sample) override;
  void AddFreespaceSample(const ion::math::Point2f& sample) override;
  std::unique_ptr<ImplicitSilhouette> Resolve() const override;

 private:
  // Quantizes the |sample| into an integer index into the solid_ and
  // freespace_ arrays.
  int Quantize(const ion::math::Point2f& sample) const;

  // Inverse of Quantize().
  ion::math::Point2f Dequantize(int index) const;

  // The size of the solid_ and freespace_ arrays.
  const ion::math::Vector2i resolution_;

  // A logical 2D array, in scanline order, indicating the presence of a solid
  // sample in a particular cell.
  std::vector<bool> solid_;

  // A logical 2D array, in scanline order, indicating the presence of a
  // freespace sample in a particular cell.
  std::vector<bool> freespace_;
};

}  // namespace baker
}  // namespace seurat

#endif  // VR_SEURAT_BAKER_FRAMEWORK_SILHOUETTE_H_
