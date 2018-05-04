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

#include "seurat/ingest/json_converter.h"

#include "gtest/gtest.h"
#include "json/json.h"
#include "seurat/api/api.pb.h"
#include "seurat/api/camera.pb.h"
#include "seurat/api/image.pb.h"
#include "seurat/api/math.pb.h"
#include "seurat/base/file_system.h"
#include "seurat/base/math_util.h"
#include "seurat/ingest/json_utils.h"
#include "seurat/ingest/json_validator.h"
#include "seurat/testing/test_flags.h"

namespace seurat {
namespace ingest {

TEST(JsonConverter, ConvertsCorrectly) {
  std::shared_ptr<base::FileSystem> file_system =
      std::make_shared<base::FileSystem>(testing::GetTestSrcdir());
  const char kJsonFilename[] =
      "com_google_seurat/seurat/ingest/testdata/manifest_v3.json";
  Json::Value root;
  EXPECT_TRUE(
      JsonUtils::ReadJsonFromFile(kJsonFilename, file_system, &root).ok());
  EXPECT_TRUE(JsonValidator::ValidateCapture(root).ok());

  api::proto::Capture capture;
  JsonConverter::ConvertCapture(root, &capture);
  EXPECT_TRUE(capture.has_headbox_center());
  EXPECT_EQ(1.0f, capture.headbox_center().x());
  EXPECT_EQ(2.0f, capture.headbox_center().y());
  EXPECT_EQ(3.0f, capture.headbox_center().z());

  EXPECT_EQ(1, capture.view_groups_size());
  api::proto::ViewGroup view_group = capture.view_groups(0);
  EXPECT_EQ(1, view_group.views_size());
  api::proto::View view = view_group.views(0);
  EXPECT_TRUE(view.has_camera());
  api::proto::Camera camera = view.camera();

  // Projective camera.
  EXPECT_TRUE(camera.has_projective());
  api::proto::ProjectiveCamera projective_camera = camera.projective();
  EXPECT_TRUE(projective_camera.has_image_size());
  EXPECT_TRUE(projective_camera.image_size().has_x());
  EXPECT_TRUE(projective_camera.image_size().has_y());
  EXPECT_EQ(450, projective_camera.image_size().x());
  EXPECT_EQ(450, projective_camera.image_size().y());
  EXPECT_TRUE(projective_camera.has_clip_from_eye());
  EXPECT_TRUE(projective_camera.has_world_from_eye());
  EXPECT_EQ(ion::math::Matrix4f::Identity(),
            base::Matrix4fFromProto(projective_camera.clip_from_eye()));
  EXPECT_EQ(2.0f * ion::math::Matrix4f::Identity(),
            base::Matrix4fFromProto(projective_camera.world_from_eye()));
  EXPECT_EQ(api::proto::DEPTH_TYPE_EYE_Z, projective_camera.depth_type());

  EXPECT_TRUE(view.has_ldi());
  api::proto::Ldi ldi = view.ldi();
  EXPECT_TRUE(ldi.has_depth_image_file());
  api::proto::DepthImageFile depth_image_file = ldi.depth_image_file();

  // Color image.
  EXPECT_TRUE(depth_image_file.has_color());
  api::proto::Image4File image4_file = depth_image_file.color();
  EXPECT_TRUE(image4_file.has_path());
  EXPECT_EQ("color.exr", image4_file.path());
  EXPECT_TRUE(image4_file.has_channel_0());
  EXPECT_EQ("R", image4_file.channel_0());
  EXPECT_TRUE(image4_file.has_channel_1());
  EXPECT_EQ("G", image4_file.channel_1());
  EXPECT_TRUE(image4_file.has_channel_2());
  EXPECT_EQ("B", image4_file.channel_2());
  EXPECT_TRUE(image4_file.has_channel_alpha());
  EXPECT_EQ("CONSTANT_ONE", image4_file.channel_alpha());

  // Depth image
  EXPECT_TRUE(depth_image_file.has_depth());
  api::proto::Image1File image1_file = depth_image_file.depth();
  EXPECT_TRUE(image1_file.has_path());
  EXPECT_EQ("depth.exr", image1_file.path());
  EXPECT_TRUE(image1_file.has_channel_0());
  EXPECT_EQ("A", image1_file.channel_0());
}

}  // namespace ingest
}  // namespace seurat
