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

#include "seurat/baker/framework/silhouette.h"

#include "ion/math/vectorutils.h"
#include "seurat/base/util.h"

namespace seurat {
namespace baker {

using ion::math::Point2f;
using ion::math::Point2i;
using ion::math::Point2ui8;
using ion::math::Vector2f;
using ion::math::Vector2ui8;

ImplicitSilhouette::ImplicitSilhouette(
    std::vector<ion::math::Point2f> solid,
    std::vector<ion::math::Point2f> freespace)
    : solid_(std::move(solid)),
      freespace_(std::move(freespace)),
      solid_tree_(solid_),
      freespace_tree_(freespace_) {}

bool ImplicitSilhouette::IsSolidAtPoint(const Point2f& point) const {
  int query_result;
  if (!solid_tree_.NNSearch(point, &query_result)) {
    // If there are no solid samples, err on the side of a completely empty
    // (freespace) region.
    return false;
  }
  const Point2f& nearest_solid = solid_[query_result];
  float distance_to_solid = ion::math::Length(nearest_solid - point);

  if (!freespace_tree_.NNSearch(point, &query_result)) {
    // If there are solid samples, but no freespace samples, the region is
    // completely solid, without any silhouette.
    return true;
  }
  const Point2f& nearest_freespace = freespace_[query_result];
  float distance_to_freespace = ion::math::Length(nearest_freespace - point);

  return distance_to_solid < distance_to_freespace;
}

namespace {

// Returns true if |p| is in [0.0, 1.0)^2.
bool InUnitSquare(const Point2f& p) {
  return p[0] < 1.0f && p[0] >= 0.0f && p[1] < 1.0f && p[1] >= 0.0f;
}

}  // namespace

void CompactSilhouetteBuffer::AddSolidSample(const Point2f& sample) {
  if (!InUnitSquare(sample)) {
    return;
  }
  solid_[Quantize(sample)] = true;
}

void CompactSilhouetteBuffer::AddFreespaceSample(const Point2f& sample) {
  if (!InUnitSquare(sample)) {
    return;
  }
  freespace_[Quantize(sample)] = true;
}

int CompactSilhouetteBuffer::Quantize(const Point2f& sample) const {
  Point2f index_f = sample * (Point2f::Zero() + Vector2f(resolution_));
  // Truncate.
  Point2i index_i = Point2i(index_f);
  return index_i[1] * resolution_[0] + index_i[0];
}

Point2f CompactSilhouetteBuffer::Dequantize(int index) const {
  Point2i index_i(index % resolution_[0], index / resolution_[0]);
  return (Point2f(index_i) + Vector2f(0.5f, 0.5f)) /
         (Point2f::Zero() + Vector2f(resolution_));
}

std::unique_ptr<ImplicitSilhouette> CompactSilhouetteBuffer::Resolve() const {
  std::vector<Point2f> solid_points;
  std::vector<Point2f> freespace_points;
  for (int i = 0; i < solid_.size(); ++i) {
    if (solid_[i] && !freespace_[i]) {
      solid_points.push_back(Dequantize(i));
    } else if (!solid_[i] && freespace_[i]) {
      freespace_points.push_back(Dequantize(i));
    }
  }
  return std::unique_ptr<ImplicitSilhouette>(new ImplicitSilhouette(
      std::move(solid_points), std::move(freespace_points)));
}

}  // namespace baker
}  // namespace seurat
