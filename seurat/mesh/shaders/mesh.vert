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

// Is this the pass for the right eye?
uniform int uIsRightEye;

// The left and right clip-from-mesh-space matrices.
uniform mat4 uClipFromMeshMatrix[2];

in highp vec3 aPosition;
in highp vec3 aTexCoord;

out highp vec3 vTexCoord;

void main() {
  vTexCoord = aTexCoord;
  gl_Position = uClipFromMeshMatrix[uIsRightEye] * vec4(aPosition, 1.0);
}
