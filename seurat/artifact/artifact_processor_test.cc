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

#include "seurat/artifact/artifact_processor.h"

#include "gtest/gtest.h"
#include "seurat/artifact/artifact.h"
#include "seurat/base/color.h"
#include "seurat/image/image.h"

namespace seurat {
namespace artifact {
namespace {

using base::Color4f;
using image::Image4f;

// Makes an artifact that has a texture with size 1x1 and the specified color
// value. It is used like a float value in
// the tests below.
Artifact MakeColorArtifact(const Color4f& color) {
  Artifact artifact;
  artifact.texture = std::make_shared<Image4f>(1, 1, color);
  return artifact;
}

// A processor that expects a ColorArtifact and adds a color value to it.
class AddProcessor : public ArtifactProcessor {
 public:
  explicit AddProcessor(const Color4f& value) : value_(value) {}
  base::Status Process(Artifact* artifact) const override {
    if (!artifact->texture) {
      return base::InvalidArgumentError(
          "Input artifact cannot be converted to texture.");
    }
    const Image4f& input_texture = *artifact->texture;
    Color4f input_color = input_texture.At(0, 0);
    artifact->texture = std::make_shared<Image4f>(1, 1, input_color + value_);
    return base::OkStatus();
  }

 private:
  const Color4f value_;
};

// A processor that expects an artifact with a texture as input and validates
// that it has a particular value.
class ValidatingProcessor : public ArtifactProcessor {
 public:
  explicit ValidatingProcessor(const Color4f& expected_color)
      : expected_color_(expected_color) {}

  base::Status Process(Artifact* artifact) const override {
    if (!artifact->texture) {
      return base::InvalidArgumentError(
          "Input artifact cannot be converted to texture.");
    }
    const Image4f& input_texture = *artifact->texture;
    Color4f actual_color = input_texture.At(0, 0);
    if (actual_color != expected_color_) {
      return base::FailedPreconditionError("Colors don't match.");
    }
    return base::OkStatus();
  }

 private:
  const Color4f expected_color_;
};

TEST(ArtifactProcessor, Sequence) {
  // Create and test the following pipeline:
  //
  // Add(1)-V(1)-Add(1)-V(2)-Add(1)-V(3)

  const Color4f kOne(1.0f, 1.0f, 1.0f, 1.0f);
  const Color4f kTwo(2.0f, 2.0f, 2.0f, 2.0f);
  const Color4f kThree(3.0f, 3.0f, 3.0f, 3.0f);

  auto sequence = std::make_shared<const ArtifactProcessorSequence>(
      std::vector<std::shared_ptr<const ArtifactProcessor>>{
          {std::make_shared<AddProcessor>(kOne),
           std::make_shared<ValidatingProcessor>(kOne),
           std::make_shared<AddProcessor>(kOne),
           std::make_shared<ValidatingProcessor>(kTwo),
           std::make_shared<AddProcessor>(kOne),
           std::make_shared<ValidatingProcessor>(kThree)}});

  Artifact artifact = MakeColorArtifact(Color4f::Zero());
  EXPECT_TRUE(sequence->Process(&artifact).ok());

  EXPECT_TRUE(artifact.texture);
  const Image4f& result_texture = *artifact.texture;
  EXPECT_EQ(Color4f(3.0f, 3.0f, 3.0f, 3.0f), result_texture.At(0, 0));
}

TEST(ArtifactProcessor, Group) {
  // Create and test the following pipeline:
  //
  //               Add(1)-V(1)
  //              /
  //         V(0)-
  //        /     \
  //       /       Add(1)-V(1)
  //      /
  //     / Add(2)-V(2)
  //     \/
  //      \
  //       Add(3)-V(3)

  const Color4f kZero(0.0f, 0.0f, 0.0f, 0.0f);
  const Color4f kOne(1.0f, 1.0f, 1.0f, 1.0f);
  const Color4f kTwo(2.0f, 2.0f, 2.0f, 2.0f);
  const Color4f kThree(3.0f, 3.0f, 3.0f, 3.0f);

  auto sequence_a = std::make_shared<const ArtifactProcessorSequence>(
      std::vector<std::shared_ptr<const ArtifactProcessor>>{
          {std::make_shared<AddProcessor>(kOne),
           std::make_shared<ValidatingProcessor>(kOne)}});

  auto group_a = std::make_shared<const ArtifactProcessorGroup>(
      std::vector<std::shared_ptr<const ArtifactProcessor>>{
          {sequence_a, sequence_a}});

  auto sequence_b = std::make_shared<const ArtifactProcessorSequence>(
      std::vector<std::shared_ptr<const ArtifactProcessor>>{
          {std::make_shared<ValidatingProcessor>(kZero), std::move(group_a)}});

  auto sequence_c = std::make_shared<const ArtifactProcessorSequence>(
      std::vector<std::shared_ptr<const ArtifactProcessor>>{
          {std::make_shared<AddProcessor>(kTwo),
           std::make_shared<ValidatingProcessor>(kTwo)}});

  auto sequence_d = std::make_shared<const ArtifactProcessorSequence>(
      std::vector<std::shared_ptr<const ArtifactProcessor>>{
          {std::make_shared<AddProcessor>(kThree),
           std::make_shared<ValidatingProcessor>(kThree)}});

  auto group_b = std::make_shared<const ArtifactProcessorGroup>(
      std::vector<std::shared_ptr<const ArtifactProcessor>>{
          {sequence_c, sequence_d}});

  auto root = std::make_shared<const ArtifactProcessorGroup>(
      std::vector<std::shared_ptr<const ArtifactProcessor>>{
          {sequence_b, group_b}});

  Artifact artifact = MakeColorArtifact(Color4f::Zero());
  EXPECT_TRUE(root->Process(&artifact).ok());

  EXPECT_TRUE(artifact.texture);
  EXPECT_EQ(Color4f::Zero(), artifact.texture->At(0, 0));
}

}  // namespace
}  // namespace artifact
}  // namespace seurat
