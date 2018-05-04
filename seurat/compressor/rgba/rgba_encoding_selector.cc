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

#include "seurat/compressor/rgba/rgba_encoding_selector.h"

#include <utility>

#include "ion/math/range.h"
#include "seurat/base/util.h"
#include "seurat/image/atlaser.h"

namespace seurat {
namespace compressor {

using image::Atlaser;
using ion::math::Point2i;
using ion::math::Range1f;
using ion::math::Vector2i;

namespace {

// Acceptable relative error in rate, distortion and Lagrange multiplier
// (lambda) computations.
constexpr float kTolerance = 0.001f;

// Maximum number of iterations.
constexpr int kMaxIterationCount = 50;

// Indicates whether lambda's value is too low or too high.
enum class LambdaState { kLow, kHigh };

struct OptimizationState {
  // Indicates whether lambda's value is too low or too high.
  LambdaState lambda_state;
  // True if all constraints are satisfied.
  bool feasible;
};

// Type of a function taking bitrate, distortion, texture encodings and selected
// encodings indices as arguments and returning the corresponding optimization
// state.
using OptimizationStateEvaluator = std::function<OptimizationState(
    float rate, float distortion,
    absl::Span<const std::vector<EncodedTexture>> encodings,
    absl::Span<const int> selection)>;

// Returns the |rate| and |distortion| values, as well as the array indices of
// the |selection| that corresponds to a given |lambda| value. The two
// dimensional array |encodings| contains the complete set of encodings for all
// textures.
void SelectEncodings(float lambda,
                     absl::Span<const std::vector<EncodedTexture>> encodings,
                     float* rate, float* distortion,
                     absl::Span<int> selection) {
  *distortion = 0.0f;
  *rate = 0.0f;
  // For each texture select the encoding that minimizes L = D + lambda R.
  for (size_t i = 0; i < encodings.size(); ++i) {
    float min_lagrangian = std::numeric_limits<float>::max();

    // Index of the selected encoding for the i-th texture.
    int selected = -1;
    for (size_t j = 0; j < encodings[i].size(); ++j) {
      float lagrangian =
          encodings[i][j].distortion + lambda * encodings[i][j].rate;
      if (lagrangian < min_lagrangian) {
        selected = j;
        min_lagrangian = lagrangian;
      }
    }
    CHECK_LE(0, selected);
    *distortion += encodings[i][selected].distortion;
    *rate += encodings[i][selected].rate;
    selection[i] = selected;
  }
}

// Returns an upper bound for the range of lambda values. The optimal lambda
// value lays between zero and this upper bound. The two dimensional array
// |encodings| contains the complete set of encodings for all textures. The
// |lambda_state_evaluator| function specifies the rate or distortion
// constraint. Also returned are the |optimization_state| and the optimal
// selection of encodings corresponding to the upper bound.
float BracketLambda(
    absl::Span<const std::vector<EncodedTexture>> encodings,
    const OptimizationStateEvaluator& optimization_state_evaluator,
    absl::Span<int> best_encodings) {
  int num_textures = encodings.size();
  float max_lambda = kTolerance;
  OptimizationState optimization_state;
  bool found_feasible = false;
  for (int i = 0; i < kMaxIterationCount; ++i) {
    float rate, distortion;
    std::vector<int> selection(num_textures);
    SelectEncodings(max_lambda, encodings, &rate, &distortion,
                    absl::MakeSpan(selection));
    optimization_state =
        optimization_state_evaluator(rate, distortion, encodings, selection);
    if (optimization_state.feasible) {
      for (int j = 0; j < num_textures; ++j) {
        best_encodings[j] = selection[j];
      }
      found_feasible = true;
    }
    if (LambdaState::kHigh == optimization_state.lambda_state) {
      break;
    }
    max_lambda *= 2.0f;
  }
  if ((LambdaState::kHigh != optimization_state.lambda_state) ||
      !found_feasible) {
    LOG(ERROR) << "Cannot satisfy the constraint imposed on content adaptive"
                  " texture resolution optimization.";
  }
  return max_lambda;
}

// Returns via |encoded_textures| the optimal selection of texture encodings,
// given the optimization objective and constraints implemented in
// |optimization_state_evaluator|. |encodings| contains the complete set of
// encodings for all textures.
void Optimize(absl::Span<const std::vector<EncodedTexture>> encodings,
              const OptimizationStateEvaluator& optimization_state_evaluator,
              absl::Span<EncodedTexture> encoded_textures) {
  int num_textures = encodings.size();
  std::vector<int> best_encodings(num_textures);
  float lower_bound = 0.0f;
  float upper_bound = BracketLambda(encodings, optimization_state_evaluator,
                                    absl::MakeSpan(best_encodings));
  // Solve for lambda using the bisection method.
  for (int iteration = 0; iteration < kMaxIterationCount; ++iteration) {
    float lambda = 0.5f * (lower_bound + upper_bound);
    if ((upper_bound - lower_bound) < kTolerance * lambda) {
      break;
    }
    float rate, distortion;
    std::vector<int> selection(encodings.size());
    SelectEncodings(lambda, encodings, &rate, &distortion,
                    absl::MakeSpan(selection));
    OptimizationState optimization_state =
        optimization_state_evaluator(rate, distortion, encodings, selection);
    if (optimization_state.feasible) best_encodings = std::move(selection);
    if (LambdaState::kLow == optimization_state.lambda_state) {
      lower_bound = lambda;
    } else if (LambdaState::kHigh == optimization_state.lambda_state) {
      upper_bound = lambda;
    }
  }
  for (int i = 0; i < num_textures; ++i) {
    encoded_textures[i] = encodings[i][best_encodings[i]];
  }
}

}  // namespace

AtlasSizeTargetSelector::AtlasSizeTargetSelector(
    std::shared_ptr<Atlaser> atlaser)
    : atlaser_(std::move(atlaser)) {}

void AtlasSizeTargetSelector::Select(
    absl::Span<const std::vector<EncodedTexture>> encodings,
    absl::Span<EncodedTexture> encoded_textures) {
  auto optimization_state_evaluator =
      [this](float rate, float distortion,
             absl::Span<const std::vector<EncodedTexture>> encodings,
             absl::Span<const int> best_encodings) {
        constexpr float kBitsPerPixel = 32.0f;
        Vector2i atlas_size_target = atlaser_->GetAtlasSizeTarget();
        float rate_target =
            kBitsPerPixel * atlas_size_target[0] * atlas_size_target[1];
        // Avoid doing layout computations if the bitrate is larger than
        // |rate_target| since in that case the constraint on the atlas size
        // cannot be satisfied.
        if (rate_target < rate) {
          return OptimizationState{LambdaState::kLow, false};
        }
        // Verify that the size of an atlas built with the selected textures
        // will not exceed |max_atlas_size|.
        int num_textures = encodings.size();
        std::vector<Vector2i> texture_sizes(num_textures);
        for (int i = 0; i < num_textures; ++i) {
          texture_sizes[i] = encodings[i][best_encodings[i]].image.GetSize();
        }
        Vector2i atlas_size;
        std::vector<Point2i> layout(num_textures);
        atlaser_->LayoutTiles(texture_sizes, &atlas_size,
                              absl::MakeSpan(layout));
        if (atlas_size[0] > atlas_size_target[0] ||
            atlas_size[1] > atlas_size_target[1]) {
          return OptimizationState{LambdaState::kLow, false};
        }
        return OptimizationState{LambdaState::kHigh, true};
      };
  Optimize(encodings, optimization_state_evaluator, encoded_textures);
}

}  // namespace compressor
}  // namespace seurat
