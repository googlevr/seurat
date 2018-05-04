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

#ifndef VR_SEURAT_IMAGE_LOW_DISCREPANCY_SAMPLING_H_
#define VR_SEURAT_IMAGE_LOW_DISCREPANCY_SAMPLING_H_

#include "ion/math/vector.h"

namespace seurat {
namespace image {

namespace low_discrepancy_sampling_internal {

float Equi(uint32 i, uint32 total) {
  return static_cast<float>(i) / static_cast<float>(total);
}

// Computes the Larcher-Pillichshammer radical inverse function introduced in:
// [Larcher and Pillichshammer] "Walsh Series Analysis of the L2-Discrepancy of
// Symmetrisized Point Sets", 2001
// The implementation is following this paper:
// [Kollig and Keller] "Efficient Multidimensional Sampling", 2002
float LarcherPillichshammer(uint32 i, uint32 scrambler = 0) {
  uint32 result = scrambler;
  uint32 v = 1 << 31;
  while (i) {
    if (i & 1) result ^= v;
    i >>= 1;
    v |= v >> 1;
  }
  return static_cast<float>(result) / static_cast<float>(0x100000000LL);
}

}  // namespace low_discrepancy_sampling_internal

// Computes the |i|-th point in the Larcher-Pillichshammer point set
// with the given |total| number of points. The point set is
// randomized using random digit scrambling with the given
// |scrambler|. A value of zero for the |scrambler| generates the
// un-scrambled point set. Larcher-Pillichshammer point sets have very
// low discrepancy and are relatively fast to compute.
ion::math::Point2f LarcherPillichshammer2D(uint32 i, uint32 total,
                                           uint32 scrambler) {
  return ion::math::Point2f(
      low_discrepancy_sampling_internal::Equi(i, total),
      low_discrepancy_sampling_internal::LarcherPillichshammer(i, scrambler));
}

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_LOW_DISCREPANCY_SAMPLING_H_
