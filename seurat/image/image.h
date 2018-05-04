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

#ifndef VR_SEURAT_IMAGE_IMAGE_H_
#define VR_SEURAT_IMAGE_IMAGE_H_

#include "seurat/base/array2d.h"
#include "seurat/base/array2d_view.h"
#include "seurat/base/color.h"

namespace seurat {
namespace image {

template <typename T>
using Image = base::Array2D<T>;

using Image1i8 = Image<base::Color1i8>;
using Image1ui8 = Image<base::Color1ui8>;
using Image1i16 = Image<base::Color1i16>;
using Image1ui16 = Image<base::Color1ui16>;
using Image1i = Image<base::Color1i>;
using Image1ui = Image<base::Color1ui>;
using Image1f = Image<base::Color1f>;
using Image1d = Image<base::Color1d>;
using Image2i8 = Image<base::Color2i8>;
using Image2ui8 = Image<base::Color2ui8>;
using Image2i16 = Image<base::Color2i16>;
using Image2ui16 = Image<base::Color2ui16>;
using Image2i = Image<base::Color2i>;
using Image2ui = Image<base::Color2ui>;
using Image2f = Image<base::Color2f>;
using Image2d = Image<base::Color2d>;
using Image3i8 = Image<base::Color3i8>;
using Image3ui8 = Image<base::Color3ui8>;
using Image3i16 = Image<base::Color3i16>;
using Image3ui16 = Image<base::Color3ui16>;
using Image3i = Image<base::Color3i>;
using Image3ui = Image<base::Color3ui>;
using Image3f = Image<base::Color3f>;
using Image3d = Image<base::Color3d>;
using Image4i8 = Image<base::Color4i8>;
using Image4ui8 = Image<base::Color4ui8>;
using Image4i16 = Image<base::Color4i16>;
using Image4ui16 = Image<base::Color4ui16>;
using Image4i = Image<base::Color4i>;
using Image4ui = Image<base::Color4ui>;
using Image4f = Image<base::Color4f>;
using Image4d = Image<base::Color4d>;

template <typename T>
using ImageView = base::Array2DView<T>;

using ImageView1i8 = ImageView<base::Color1i8>;
using ImageView1ui8 = ImageView<base::Color1ui8>;
using ImageView1i16 = ImageView<base::Color1i16>;
using ImageView1ui16 = ImageView<base::Color1ui16>;
using ImageView1i = ImageView<base::Color1i>;
using ImageView1ui = ImageView<base::Color1ui>;
using ImageView1f = ImageView<base::Color1f>;
using ImageView1d = ImageView<base::Color1d>;
using ImageView2i8 = ImageView<base::Color2i8>;
using ImageView2ui8 = ImageView<base::Color2ui8>;
using ImageView2i16 = ImageView<base::Color2i16>;
using ImageView2ui16 = ImageView<base::Color2ui16>;
using ImageView2i = ImageView<base::Color2i>;
using ImageView2ui = ImageView<base::Color2ui>;
using ImageView2f = ImageView<base::Color2f>;
using ImageView2d = ImageView<base::Color2d>;
using ImageView3i8 = ImageView<base::Color3i8>;
using ImageView3ui8 = ImageView<base::Color3ui8>;
using ImageView3i16 = ImageView<base::Color3i16>;
using ImageView3ui16 = ImageView<base::Color3ui16>;
using ImageView3i = ImageView<base::Color3i>;
using ImageView3ui = ImageView<base::Color3ui>;
using ImageView3f = ImageView<base::Color3f>;
using ImageView3d = ImageView<base::Color3d>;
using ImageView4i8 = ImageView<base::Color4i8>;
using ImageView4ui8 = ImageView<base::Color4ui8>;
using ImageView4i16 = ImageView<base::Color4i16>;
using ImageView4ui16 = ImageView<base::Color4ui16>;
using ImageView4i = ImageView<base::Color4i>;
using ImageView4ui = ImageView<base::Color4ui>;
using ImageView4f = ImageView<base::Color4f>;
using ImageView4d = ImageView<base::Color4d>;

template <typename T>
using MutableImageView = base::MutableArray2DView<T>;

using MutableImageView1i8 = MutableImageView<base::Color1i8>;
using MutableImageView1ui8 = MutableImageView<base::Color1ui8>;
using MutableImageView1i16 = MutableImageView<base::Color1i16>;
using MutableImageView1ui16 = MutableImageView<base::Color1ui16>;
using MutableImageView1i = MutableImageView<base::Color1i>;
using MutableImageView1ui = MutableImageView<base::Color1ui>;
using MutableImageView1f = MutableImageView<base::Color1f>;
using MutableImageView1d = MutableImageView<base::Color1d>;
using MutableImageView2i8 = MutableImageView<base::Color2i8>;
using MutableImageView2ui8 = MutableImageView<base::Color2ui8>;
using MutableImageView2i16 = MutableImageView<base::Color2i16>;
using MutableImageView2ui16 = MutableImageView<base::Color2ui16>;
using MutableImageView2i = MutableImageView<base::Color2i>;
using MutableImageView2ui = MutableImageView<base::Color2ui>;
using MutableImageView2f = MutableImageView<base::Color2f>;
using MutableImageView2d = MutableImageView<base::Color2d>;
using MutableImageView3i8 = MutableImageView<base::Color3i8>;
using MutableImageView3ui8 = MutableImageView<base::Color3ui8>;
using MutableImageView3i16 = MutableImageView<base::Color3i16>;
using MutableImageView3ui16 = MutableImageView<base::Color3ui16>;
using MutableImageView3i = MutableImageView<base::Color3i>;
using MutableImageView3ui = MutableImageView<base::Color3ui>;
using MutableImageView3f = MutableImageView<base::Color3f>;
using MutableImageView3d = MutableImageView<base::Color3d>;
using MutableImageView4i8 = MutableImageView<base::Color4i8>;
using MutableImageView4ui8 = MutableImageView<base::Color4ui8>;
using MutableImageView4i16 = MutableImageView<base::Color4i16>;
using MutableImageView4ui16 = MutableImageView<base::Color4ui16>;
using MutableImageView4i = MutableImageView<base::Color4i>;
using MutableImageView4ui = MutableImageView<base::Color4ui>;
using MutableImageView4f = MutableImageView<base::Color4f>;
using MutableImageView4d = MutableImageView<base::Color4d>;

}  // namespace image
}  // namespace seurat

#endif  // VR_SEURAT_IMAGE_IMAGE_H_
