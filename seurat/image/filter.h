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

#ifndef VR_SEURAT_IMAGE_FILTER_H_
#define VR_SEURAT_IMAGE_FILTER_H_

#include "ion/math/range.h"
#include "ion/math/vector.h"
#include "ion/math/vectorutils.h"

namespace seurat {
namespace image {

class Filter {
 public:
  // Make destructor virtual.
  virtual ~Filter() = default;

  // Evaluate the filter for the given position relative to it's center.
  virtual float Eval(const ion::math::Point2f& position) const = 0;

  // Return the radius of the filter. This is a bounding circle around the
  // support of the filter.
  virtual float GetRadius() const = 0;

  // Returns the bounding box of the filter support. Coordinates are relative to
  // the center of the filter.
  virtual ion::math::Range2f GetRange() const = 0;
};

class BoxFilter : public Filter {
 public:
  // Constructs a box filter with the given |size|.
  explicit BoxFilter(const ion::math::Vector2f& size) {
    const ion::math::Vector2f half_size = size * 0.5f;
    radius_ = ion::math::Length(half_size);
    range_ = ion::math::Range2f(ion::math::Point2f::Zero() - half_size,
                                ion::math::Point2f::Zero() + half_size);
  }
  // Constructs a box filter with a unit square as the support.
  BoxFilter() : BoxFilter(ion::math::Vector2f(1.0f, 1.0f)) {}

  float Eval(const ion::math::Point2f& position) const override {
    return range_.ContainsPoint(position);
  }

  float GetRadius() const override { return radius_; }

  ion::math::Range2f GetRange() const override { return range_; }

 private:
  ion::math::Range2f range_;
  float radius_;
};

class BSplineFilter : public Filter {
 public:
  float Eval(const ion::math::Point2f& position) const override {
    const float d = ion::math::Length(position - ion::math::Point2f::Zero());
    if (d > 2.0f) {
      return 0.0f;
    } else if (d < 1.0f) {
      const float t = 1.0f - d;
      return ((((-3.0f * t) + 3.0f) * t + 3.0f) * t + 1.0f) / 6.0f;
    } else {
      const float t = 2.0f - d;
      return t * t * t / 6.0f;
    }
  }

  float GetRadius() const override { return 2.0f; }

  ion::math::Range2f GetRange() const override {
    return ion::math::Range2f(ion::math::Point2f(-2.0f, -2.0f),
                              ion::math::Point2f(2.0f, 2.0f));
  }
};

class MitchellFilter : public Filter {
 public:
  // Constructs a Mitchell filter with the given parameters |B| and |C|.
  // B = 0, C = 1 is the cubic B-spline.
  // B = 0, is the family of cardinal splines.
  // B = 0, C = 0.5 is the Catmull-Rom spline.
  // The authors of the original papaer suggest B + 2C = 1 as good parameters.
  // In particular B = C = 1/3.
  MitchellFilter(float b, float c) : b_(b), c_(c) {}

  // Default constructor, B = C = 1/3.
  MitchellFilter() : b_(1.0f / 3.0f), c_(1.0f / 3.0f) {}

  float Eval(const ion::math::Point2f& position) const override {
    const float d = ion::math::Length(position - ion::math::Point2f::Zero());
    const float d2 = d * d;
    const float d3 = d2 * d;
    if (d < 1.0f) {
      return ((12.0f - 9.0f * b_ - 6.0f * c_) * d3 +    //
              (-18.0f + 12.0f * b_ + 6.0f * c_) * d2 +  //
              (6.0f - 2.0f * b_)) /                     //
             6.0f;
    } else if ((d >= 1.0f) && (d < 2.0f)) {
      return ((-b_ - 6.0f * c_) * d3 +          //
              (6.0f * b_ + 30.0f * c_) * d2 +   //
              (-12.0f * b_ - 48.0f * c_) * d +  //
              (8.0f * b_ + 24.0f * c_)) /       //
             6.0f;
    } else {
      return 0.0f;
    }
  }

  float GetRadius() const override { return 2.0f; }

  ion::math::Range2f GetRange() const override {
    return ion::math::Range2f(ion::math::Point2f(-2.0f, -2.0f),
                              ion::math::Point2f(2.0f, 2.0f));
  }

 private:
  // The two constants B and C defining the filter's shape.
  const float b_;
  const float c_;
};

class GaussianFilter : public Filter {
 public:
  // Constructs a gaussian filter with the given standard deviation |sigma| and
  // |radius|. The gaussian is truncated at distance |radius| and shifted such
  // that f(radius) = 0. As a rule of thumb, radius = 3 * sigma is a reasonable
  // choice for the cut-off.
  explicit GaussianFilter(float sigma, float radius)
      : radius_(radius),
        a_(1.0f / (std::sqrt(2 * M_PI) * sigma)),
        b_(-1.0f / (2.0f * sigma * sigma)),
        c_(-a_ * std::exp(radius * radius * b_)) {}

  float Eval(const ion::math::Point2f& position) const override {
    const float d2 =
        ion::math::LengthSquared(position - ion::math::Point2f::Zero());
    if (d2 >= radius_ * radius_) {
      return 0.0f;
    }
    return a_ * std::exp(d2 * b_) + c_;
  }

  float GetRadius() const override { return radius_; }

  ion::math::Range2f GetRange() const override {
    return ion::math::Range2f(ion::math::Point2f(-radius_, -radius_),
                              ion::math::Point2f(radius_, radius_));
  }

 private:
  const float radius_;
  const float a_;
  const float b_;
  const float c_;
};

class Wendland31Filter : public Filter {
 public:
  // Constructs a Wendland filter with the given |radius|, which is a
  // bell-shaped filter similar to a Gaussian. The filter belongs to a doubly
  // indexed family of filters and the indices of this particular one are (3,1).
  //
  // This filter has a number of properties that are useful for antialiasing:
  // - it is radially symmetric.
  // - it is finitely supported (zero for points farther away from the origin
  // than |radius|).
  // - it is C2 continuous at r = 0 and at r = |radius| (and Cinfinity
  // everywhere else.
  // - it has a Fourier transform that is positive everywhere, i.e. similar to a
  // non-truncated Gaussian filter. B-spline basis filters, although finitely
  // supported and non-negative everywhere have Fourier transforms with negative
  // lobes.
  // - it is simpler to evaluate than either the Mitchell or cubic B-spline
  // basis filters.
  explicit Wendland31Filter(float radius)
      : radius_(radius),
        radius_squared_(radius * radius),
        inv_radius_squared_(1.0f / radius_squared_) {}

  float Eval(const ion::math::Point2f& position) const override {
    const float d2 =
        ion::math::LengthSquared(position - ion::math::Point2f::Zero());
    if (d2 >= radius_squared_) {
      return 0.0f;
    }
    float r = std::sqrt(d2 * inv_radius_squared_);
    float omr = 1.0f - r;
    float omr2 = omr * omr;
    return omr2 * omr2 * (1.0f + 4.0f * r);
  }

  float GetRadius() const override { return radius_; }

  ion::math::Range2f GetRange() const override {
    return ion::math::Range2f(ion::math::Point2f(-radius_, -radius_),
                              ion::math::Point2f(radius_, radius_));
  }

 private:
  const float radius_;
  const float radius_squared_;
  const float inv_radius_squared_;
};

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_FILTER_H_
