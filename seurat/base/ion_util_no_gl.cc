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

#include "seurat/base/ion_util_no_gl.h"

#include "ion/base/datacontainer.h"
#include "ion/gfx/image.h"
#include "ion/gfx/sampler.h"
#include "ion/gfx/shaderinputregistry.h"
#include "ion/gfx/texture.h"
#include "ion/gfxutils/shadermanager.h"
#include "ion/gfxutils/shadersourcecomposer.h"
#include "ion/math/transformutils.h"

namespace seurat {
namespace base {
namespace {

// This utility function constructs a Point3d representing one of the vertices
// of a Range3f.  The vertices are indexed arbitrarily on the range [0, 8).
template <typename T>
ion::math::Point<3, T> SelectCorner(const ion::math::Range<3, T>& box, int i) {
  return ion::math::Point<3, T>(
      (i & 0x1) == 0 ? box.GetMinPoint()[0] : box.GetMaxPoint()[0],
      (i & 0x2) == 0 ? box.GetMinPoint()[1] : box.GetMaxPoint()[1],
      (i & 0x4) == 0 ? box.GetMinPoint()[2] : box.GetMaxPoint()[2]);
}

}  // namespace

template <typename T>
ion::math::Range<3, T> ProjectAABB(const ion::math::Matrix<4, T>& transform,
                                   const ion::math::Range<3, T>& aabb) {
  ion::math::Range<3, T> new_aabb;
  for (int i = 0; i < 8; ++i) {
    new_aabb.ExtendByPoint(
        ion::math::ProjectPoint(transform, SelectCorner(aabb, i)));
  }
  return new_aabb;
}

template ion::math::Range<3, float> ProjectAABB(
    const ion::math::Matrix<4, float>& transform,
    const ion::math::Range<3, float>& aabb);
template ion::math::Range<3, double> ProjectAABB(
    const ion::math::Matrix<4, double>& transform,
    const ion::math::Range<3, double>& aabb);

ion::gfx::ImagePtr CreateImage(ion::gfx::Image::Format format,
                               const ion::math::Vector2i& dimensions) {
  ion::gfx::ImagePtr image(new ion::gfx::Image());
  ion::base::DataContainerPtr container =
      ion::base::DataContainer::CreateOverAllocated<uint8>(
          ion::gfx::Image::ComputeDataSize(format, dimensions[0],
                                           dimensions[1]),
          NULL, image->GetAllocator());
  image->Set(format, dimensions[0], dimensions[1], container);
  return image;
}

ion::gfx::TexturePtr CreateTexture(ion::gfx::ImagePtr image) {
  ion::gfx::TexturePtr texture(new ion::gfx::Texture);
  ion::gfx::SamplerPtr sampler(new ion::gfx::Sampler);
  texture->SetImage(0U, image);
  texture->SetSampler(sampler);
  return texture;
}

bool CompareImagesEqual(const ion::gfx::ImagePtr& lhs,
                        const ion::gfx::ImagePtr& rhs) {
  if (lhs.Get() == nullptr) {
    return rhs.Get() == nullptr;
  }
  if (rhs.Get() == nullptr) {
    return false;
  }
  if (lhs->GetFormat() != rhs->GetFormat()) {
    return false;
  }
  if (lhs->GetDimensions() != rhs->GetDimensions()) {
    return false;
  }
  if (lhs->GetDepth() != rhs->GetDepth()) {
    return false;
  }
  if (lhs->GetDataSize() != rhs->GetDataSize()) {
    return false;
  }
  return std::equal(lhs->GetData()->GetData<uint8>(),
                    lhs->GetData()->GetData<uint8>() + lhs->GetDataSize(),
                    rhs->GetData()->GetData<uint8>());
}

ion::gfx::ShaderProgramPtr CreateShaderProgram(
    const ion::gfxutils::ShaderManagerPtr& shader_manager,
    const ion::gfx::ShaderInputRegistryPtr& shader_registry,
    const string& shader_program_name, const string& vertex_shader_path,
    const string& fragment_shader_path) {
  // Prepend shaders compiled through CreateShaderProgram with the proper
  // #version directive.  We choose GLSL 3.30 / GLSL ES 3.00 because:
  //
  // * They both support parameter and precision qualifiers.
  // * GLSL 3.30 is supported by OpenGL 3.3, which has rough feature parity with
  //   OpenGL ES 2.0.
  //
  // This minimizes the amount of shader rewriting that has to happens between
  // OpenGL and OpenGL ES shaders.
  static const auto kPrependShaderVersion = [](const std::string& source) {
#if !defined(ION_GFX_OGLES20)
    static const std::string kPrefix = std::string(
        "// Prefixed by vr_softshader::base::CreateShaderProgram\n"
        "#version 330\n"
        "#line 1\n");
#else
    static const std::string kPrefix = std::string(
        "// Prefixed by vr_softshader::base::CreateShaderProgram\n"
        "#version 300 es\n"
        "#line 1\n");
#endif
    return kPrefix + source;
  };

  return shader_manager->CreateShaderProgram(
      shader_program_name, shader_registry,
      ion::gfxutils::ShaderSourceComposerPtr(new ion::gfxutils::FilterComposer(
          ion::gfxutils::ShaderSourceComposerPtr(
              new ion::gfxutils::ZipAssetComposer(vertex_shader_path, false)),
          kPrependShaderVersion)),
      ion::gfxutils::ShaderSourceComposerPtr(new ion::gfxutils::FilterComposer(
          ion::gfxutils::ShaderSourceComposerPtr(
              new ion::gfxutils::ZipAssetComposer(fragment_shader_path, false)),
          kPrependShaderVersion)));
}

}  // namespace base
}  // namespace seurat
