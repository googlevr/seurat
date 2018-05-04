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

#include "seurat/artifact/evaluation/render_sim.h"

#include "ion/math/matrix.h"
#include "ion/math/transformutils.h"
#include "ion/math/vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "seurat/base/color.h"
#include "seurat/base/projective_camera.h"
#include "seurat/base/projective_camera_util.h"
#include "seurat/geometry/mesh.h"
#include "seurat/image/image.h"
#include "seurat/image/image_util.h"
#include "seurat/testing/ion_test_utils.h"

namespace seurat {
namespace artifact {
namespace {

using base::Color4f;
using base::ProjectiveCamera;
using geometry::Mesh;
using image::Image4f;
using ion::math::Point3f;
using ion::math::Vector2i;
using ion::math::Vector3f;

TEST(RenderSimTest, TestBilinearFiltering) {
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kGreen(0.0f, 1.0f, 0.0f, 1.0f);
  const Color4f kBlue(0.0f, 0.0f, 1.0f, 1.0f);
  const Color4f kBlack(0.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kWhite(1.0f, 1.0f, 1.0f, 1.0f);

  const bool kEnableZBuffer = false;
  RenderSim render_sim(kEnableZBuffer);
  // Add a unit-size quad at z = -1.
  //
  // Note that the left half (small 'x') of the 4x2 texture will correspond to a
  // square around one texture sample.  So, when rendered, the corner nearest to
  // (0, 0, -1) should be solid, and the corner nearest to (1, 1, -1) should
  // have a linear gradient approaching an even mix of all 4 texels of the left
  // half of the texture.
  //
  // ^
  // |green  kBlack
  // 3----2
  // |    |
  // |    |
  // 0----1-->
  // red    blue
  Mesh mesh(1);
  mesh.AppendVertex(Point3f(0.0f, 0.0f, -1.0f), {{Point3f(0.0f, 0.0f, 1.0f)}});
  mesh.AppendVertex(Point3f(1.0f, 0.0f, -1.0f), {{Point3f(0.25f, 0.0f, 1.0f)}});
  mesh.AppendVertex(Point3f(1.0f, 1.0f, -1.0f), {{Point3f(0.25f, 0.5f, 1.0f)}});
  mesh.AppendVertex(Point3f(0.0f, 1.0f, -1.0f), {{Point3f(0.0f, 0.5f, 1.0f)}});
  mesh.AppendTriangle({{0, 1, 2}});
  mesh.AppendTriangle({{0, 2, 3}});

  Image4f texture(4, 2);
  texture.Fill(kWhite);
  texture.At(0, 0) = kRed;
  texture.At(1, 0) = kBlue;
  texture.At(0, 1) = kGreen;
  texture.At(1, 1) = kBlack;

  render_sim.Build(mesh, &texture);

  EXPECT_EQ(Color4f::Zero(), render_sim.TraceRay(Point3f(0.5f, 0.5f, -1.1f),
                                                 Vector3f(0.0f, 0.0f, -1.0f)));

  // The subregion of the square with (x, y) in [0, 0.5)^2 should be red.
  EXPECT_VECTOR_NEAR(kRed,
                     render_sim.TraceRay(Point3f(0.01f, 0.01f, 0.0f),
                                         Vector3f(0.0f, 0.0f, -0.1f)),
                     1e-5f);
  EXPECT_VECTOR_NEAR(kRed,
                     render_sim.TraceRay(Point3f(0.25f, 0.25f, 0.0f),
                                         Vector3f(0.0f, 0.0f, -0.1f)),
                     1e-5f);
  EXPECT_VECTOR_NEAR(kRed,
                     render_sim.TraceRay(Point3f(0.01f, 0.49f, 0.0f),
                                         Vector3f(0.0f, 0.0f, -0.1f)),
                     1e-5f);
  EXPECT_VECTOR_NEAR(kRed,
                     render_sim.TraceRay(Point3f(0.01f, 0.49f, 0.0f),
                                         Vector3f(0.0f, 0.0f, -0.1f)),
                     1e-5f);

  // A point along the bottom edge.  It should be 3/4 red and 1/4 blue.
  EXPECT_VECTOR_NEAR(kRed * 3.0f / 4.0f + kBlue * 1.0f / 4.0f,
                     render_sim.TraceRay(Point3f(0.75f, 1.0e-4f, 0.0f),
                                         Vector3f(0.0f, 0.0f, -0.1f)),
                     1e-5f);

  // A point along the left edge.  It should be 3/4 red and 1/4 green.
  EXPECT_VECTOR_NEAR(kRed * 3.0f / 4.0f + kGreen * 1.0f / 4.0f,
                     render_sim.TraceRay(Point3f(1e-4f, 0.75f, 0.0f),
                                         Vector3f(0.0f, 0.0f, -0.1f)),
                     1e-5f);
}

TEST(RenderSimTest, TestRender) {
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kGreen(0.0f, 1.0f, 0.0f, 1.0f);
  const Color4f kBlue(0.0f, 0.0f, 1.0f, 1.0f);
  const Color4f kBlack(0.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kWhite(1.0f, 1.0f, 1.0f, 1.0f);

  const int kThreadCount = 3;
  const bool kEnableZBuffer = false;
  RenderSim render_sim(kThreadCount, kEnableZBuffer);
  // Add a 2x2 quad at z = -1.
  Mesh mesh(1);
  mesh.AppendVertex(Point3f(0.0f, -1.0f, -1.0f), {{Point3f(0.0f, 0.0f, 1.0f)}});
  mesh.AppendVertex(Point3f(1.0f, -1.0f, -1.0f),
                    {{Point3f(0.25f, 0.0f, 1.0f)}});
  mesh.AppendVertex(Point3f(1.0f, 1.0f, -1.0f), {{Point3f(0.25f, 0.5f, 1.0f)}});
  mesh.AppendVertex(Point3f(-1.0f, 1.0f, -1.0f), {{Point3f(0.0f, 0.5f, 1.0f)}});
  mesh.AppendTriangle({{0, 1, 2}});
  mesh.AppendTriangle({{0, 2, 3}});

  Image4f texture(4, 2);
  texture.Fill(kWhite);
  texture.At(0, 0) = kRed;
  texture.At(1, 0) = kBlue;
  texture.At(0, 1) = kGreen;
  texture.At(1, 1) = kBlack;

  render_sim.Build(mesh, &texture);

  const Vector2i kImageSize(5, 5);
  const float kNear = 0.5f;
  const float kFar = 2.0f;
  ProjectiveCamera camera(kImageSize,
                          ion::math::PerspectiveMatrixFromFrustum(
                              -kNear, kNear, -kNear, kNear, kNear, kFar),
                          ion::math::Matrix4f::Identity());

  Image4f rendered;
  render_sim.Render(camera, &rendered);
  EXPECT_VECTOR_FLOAT_EQ(kImageSize, rendered.GetSize());

  // The center pixel should be a mix of red, green, and blue.
  EXPECT_LT(0, rendered.At(3, 3)[0]);
  EXPECT_LT(0, rendered.At(3, 3)[1]);
  EXPECT_LT(0, rendered.At(3, 3)[2]);
  EXPECT_EQ(1.0f, rendered.At(3, 3)[3]);
}

TEST(RenderSimTest, TestPrimitiveOrder) {
  const Color4f kRed(1.0f, 0.0f, 0.0f, 1.0f);
  const Color4f kGreen(0.0f, 1.0f, 0.0f, 1.0f);
  const Color4f kBlue(0.0f, 0.0f, 1.0f, 1.0f);

  // Add a sequence of unit-size quads in the following order:
  //  * red at z = -1
  //  * green at z = -3
  //  * blue at z = -2

  Image4f texture(3, 1);
  texture.At(0, 0) = kRed;
  texture.At(1, 0) = kGreen;
  texture.At(2, 0) = kBlue;

  Mesh mesh(1);

  const Point3f red_uv(1.0f / 6.0f, 0.5f, 1.0f);
  mesh.AppendVertex(Point3f(0.0f, 0.0f, -1.0f), {{red_uv}});
  mesh.AppendVertex(Point3f(1.0f, 0.0f, -1.0f), {{red_uv}});
  mesh.AppendVertex(Point3f(1.0f, 1.0f, -1.0f), {{red_uv}});
  mesh.AppendVertex(Point3f(0.0f, 1.0f, -1.0f), {{red_uv}});
  mesh.AppendTriangle({{0, 1, 2}});
  mesh.AppendTriangle({{0, 2, 3}});

  const Point3f green_uv(3.0f / 6.0f, 0.5f, 1.0f);
  mesh.AppendVertex(Point3f(0.0f, 0.0f, -3.0f), {{green_uv}});
  mesh.AppendVertex(Point3f(1.0f, 0.0f, -3.0f), {{green_uv}});
  mesh.AppendVertex(Point3f(1.0f, 1.0f, -3.0f), {{green_uv}});
  mesh.AppendVertex(Point3f(0.0f, 1.0f, -3.0f), {{green_uv}});
  mesh.AppendTriangle({{4, 5, 6}});
  mesh.AppendTriangle({{4, 6, 7}});

  const Point3f blue_uv(5.0f / 6.0f, 0.5f, 1.0f);
  mesh.AppendVertex(Point3f(0.0f, 0.0f, -2.0f), {{blue_uv}});
  mesh.AppendVertex(Point3f(1.0f, 0.0f, -2.0f), {{blue_uv}});
  mesh.AppendVertex(Point3f(1.0f, 1.0f, -2.0f), {{blue_uv}});
  mesh.AppendVertex(Point3f(0.0f, 1.0f, -2.0f), {{blue_uv}});
  mesh.AppendTriangle({{8, 9, 10}});
  mesh.AppendTriangle({{8, 10, 11}});

  const int kThreadCount = 3;

  // Trace rays which run perpendicular to all quads and verify that the
  // resulting color indicates the desired ordering.
  //
  // Note that rays are perturbed from the center to also test correct filtering
  // behavior.
  {
    // Test without a z-buffer.
    const bool kEnableZBuffer = false;
    RenderSim render_sim(kThreadCount, kEnableZBuffer);
    render_sim.Build(mesh, &texture);
    // This ray will hit all 3 squares, so it should be blue, since blue will
    // render last (even though the ray hits it second!).
    EXPECT_VECTOR_NEAR(kBlue,
                       render_sim.TraceRay(Point3f(0.1f, 0.3f, 0.0f),
                                           Vector3f(0.0f, 0.0f, -1.0f)),
                       1e-5f);
    // This ray will hit the green and blue squares, so it should be blue, since
    // blue will render last.
    EXPECT_VECTOR_NEAR(kBlue,
                       render_sim.TraceRay(Point3f(0.9f, 0.4f, -1.5f),
                                           Vector3f(0.0f, 0.0f, -1.0f)),
                       1e-5f);
    // This ray will hit the blue and red squares, so it should be blue, since
    // blue will render last.
    EXPECT_VECTOR_NEAR(kBlue,
                       render_sim.TraceRay(Point3f(0.2f, 0.8f, -2.5f),
                                           Vector3f(0.0f, 0.0f, 1.0f)),
                       1e-5f);
  }

  {
    // Test with a z-buffer.
    const bool kEnableZBuffer = true;
    RenderSim render_sim(kThreadCount, kEnableZBuffer);
    render_sim.Build(mesh, &texture);
    // This ray will hit all 3 squares, so it should be red, since red will
    // render last.
    EXPECT_VECTOR_NEAR(kRed,
                       render_sim.TraceRay(Point3f(0.1f, 0.3f, 0.0f),
                                           Vector3f(0.0f, 0.0f, -1.0f)),
                       1e-5f);
    // This ray will hit all 3 squares from the other side, so it should be
    // green, since green will render last.
    EXPECT_VECTOR_NEAR(kGreen,
                       render_sim.TraceRay(Point3f(0.9f, 0.4f, -3.5f),
                                           Vector3f(0.0f, 0.0f, 1.0f)),
                       1e-5f);
    // This ray will hit the blue and red squares, so it should be blue, since
    // blue will render last.
    EXPECT_VECTOR_NEAR(kBlue,
                       render_sim.TraceRay(Point3f(0.2f, 0.8f, -2.5f),
                                           Vector3f(0.0f, 0.0f, 1.0f)),
                       1e-5f);
  }
}

}  // namespace
}  // namespace artifact
}  // namespace seurat
