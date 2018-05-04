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

precision mediump float;

// Incoming interpolated vertex color.
in vec4 vColor;

// Output color.
out vec4 frag_color;

void main() {
  vec2 position = gl_PointCoord - vec2(0.5, 0.5);
  float sqr_distance = dot(position, position);
  // Discard fragments outside the inscribed circle of the point splat.
  if (sqr_distance > 0.25) discard;
  // Compute a weight that falls of to zero towards the edge of the splat
  // circle. The resulting images using this cheap weighting function are very
  // close to a gaussian weight.
  float weight = 1.0 - (sqr_distance * 4.0);
  frag_color = vec4(vColor * weight);
}
