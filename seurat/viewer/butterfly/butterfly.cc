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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ion/gfx/graphicsmanager.h"
#include "ion/gfx/renderer.h"
#include "ion/gfx/shaderinputregistry.h"
#include "ion/gfxutils/shadermanager.h"
#include "ion/math/matrix.h"
#include "ion/math/transformutils.h"
#include "ion/port/memorymappedfile.h"
#include "ion/remote/calltracehandler.h"
#include "ion/remote/nodegraphhandler.h"
#include "ion/remote/remoteserver.h"
#include "ion/remote/resourcehandler.h"
#include "ion/remote/settinghandler.h"
#include "ion/remote/shaderhandler.h"
#include "ion/remote/tracinghandler.h"
#include "include/GLFW/glfw3.h"
#include "absl/strings/str_cat.h"
#include "seurat/base/reporting.h"
#include "seurat/base/structured_io.h"
#include "seurat/component/component.h"
#include "seurat/component/ice_io.h"
#include "seurat/component/renderable.h"
#include "seurat/viewer/butterfly/viewer_camera.h"
#include "seurat/viewer/rendering_pipeline.h"
#include "seurat/viewer/scene.h"

using RenderingContext = seurat::component::Renderable::RenderingContext;
using ViewState = seurat::component::Renderable::ViewState;
using ion::math::Point2d;
using ion::math::Point2i;
using ion::math::Range2i;
using ion::math::Vector2d;
using ion::math::Vector3f;
using seurat::component::Component;
using seurat::viewer::ViewerCamera;

const int kWindowWidth = 1024;
const int kWindowHeight = 1024;

// User data that is attached to the GLFW window to allow stateful callbacks.
struct WindowUserData {
  // Coordinates of the curser at the previous mouse event.
  Point2d previous_cursor_position;
  // The viewer camera.
  ViewerCamera viewer_camera;
};

void ErrorCallback(int error, const char* message) {
  seurat::base::SeuratError(absl::StrCat("OpenGL error: ", message));
}

static void KeyboardCallback(GLFWwindow* window, int key, int scancode,
                             int action, int mods) {
  const float kTranslateSpeed = 0.02f;
  WindowUserData* user_data =
      static_cast<WindowUserData*>(glfwGetWindowUserPointer(window));

  std::map<int, Vector3f> direction_from_key;
  direction_from_key[GLFW_KEY_D] = Vector3f(1.0f, 0.0f, 0.0f);
  direction_from_key[GLFW_KEY_A] = Vector3f(-1.0f, 0.0f, 0.0f);
  direction_from_key[GLFW_KEY_Q] = Vector3f(0.0f, 1.0f, 0.0f);
  direction_from_key[GLFW_KEY_Z] = Vector3f(0.0f, -1.0f, 0.0f);
  direction_from_key[GLFW_KEY_S] = Vector3f(0.0f, 0.0f, 1.0f);
  direction_from_key[GLFW_KEY_W] = Vector3f(0.0f, 0.0f, -1.0f);

  // Only respond to key-press and key-repeat. Don't do anything on key release.
  if (action != GLFW_PRESS && action != GLFW_REPEAT) {
    return;
  }

  switch (key) {
    case GLFW_KEY_ESCAPE:  // Close window and exit.
      glfwSetWindowShouldClose(window, GLFW_TRUE);
      break;
    case GLFW_KEY_TAB:  // Reset camera position to origin.
      user_data->viewer_camera.Reset();
      break;
    case GLFW_KEY_D:  // Move right
    case GLFW_KEY_A:  // Move left
    case GLFW_KEY_Q:  // Move up
    case GLFW_KEY_Z:  // Move down
    case GLFW_KEY_S:  // Move backward
    case GLFW_KEY_W:  // Move forward
      user_data->viewer_camera.Translate(direction_from_key[key] *
                                         kTranslateSpeed);
      break;
  }
}

void MouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
  WindowUserData* user_data =
      static_cast<WindowUserData*>(glfwGetWindowUserPointer(window));
  if (action == GLFW_PRESS) {
    // Save the cursor position on any mouse button down event.
    glfwGetCursorPos(window, &user_data->previous_cursor_position[0],
                     &user_data->previous_cursor_position[1]);
  }
}

static void CursorPositionCallback(GLFWwindow* window, double position_x,
                                   double position_y) {
  const float kRotateSpeed = 0.25f;
  const float kTranslateSpeed = 0.01f;

  WindowUserData* user_data =
      static_cast<WindowUserData*>(glfwGetWindowUserPointer(window));
  const Point2d cursor_position(position_x, position_y);
  const Vector2d cursor_delta =
      cursor_position - user_data->previous_cursor_position;
  user_data->previous_cursor_position = cursor_position;

  if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
    // Mouse move with LMB pressed: Rotate camera.
    user_data->viewer_camera.Rotate(
        ion::math::Anglef(),
        ion::math::Anglef::FromDegrees(cursor_delta[1] * kRotateSpeed),
        ion::math::Anglef::FromDegrees(cursor_delta[0] * kRotateSpeed));
  } else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) ==
             GLFW_PRESS) {
    // Mouse move with MMB pressed: Translate camera in eye-space.
    user_data->viewer_camera.Translate(
        ion::math::Vector3f(cursor_delta[0] * kTranslateSpeed,
                            -cursor_delta[1] * kTranslateSpeed, 0.0f));
  } else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) ==
             GLFW_PRESS) {
    // Mouse move with RMB pressed: Move camera forward and backward in
    // eye-space.
    user_data->viewer_camera.Translate(
        ion::math::Vector3f(0.0f, 0.0f, cursor_delta[1] * kTranslateSpeed));
  }
}

// Computes the view state for the given |viewer_camera|.
ViewState ComputeViewState(const ViewerCamera& viewer_camera) {
  const ion::math::Matrix4f eye_from_world = viewer_camera.GetEyeFromWorld();
  const ion::math::Matrix4f clip_from_eye = viewer_camera.GetClipFromEye();
  ViewState view_state;
  view_state.left_eye_from_component_matrix = eye_from_world;
  view_state.right_eye_from_component_matrix = eye_from_world;
  view_state.left_clip_from_component_matrix = clip_from_eye * eye_from_world;
  view_state.right_clip_from_component_matrix = clip_from_eye * eye_from_world;
  return view_state;
}

// Loads and returns the component stored in the file with the given |filename|.
std::unique_ptr<const Component> LoadComponent(const std::string& filename) {
  ion::port::MemoryMappedFile file(filename);
  if (file.GetData() == nullptr) {
    seurat::base::SeuratError(
        absl::StrCat("cannot open ", filename, " for reading"));
    exit(EXIT_FAILURE);
  }
  absl::string_view file_string(static_cast<const char*>(file.GetData()),
                                file.GetLength());
  seurat::base::ArrayByteSource byte_source(file_string);
  seurat::base::StructureSource source(&byte_source);
  return seurat::component::ReadIce(&source);
}

// Creates, initializes and returns a Seurat renderer.
std::unique_ptr<seurat::viewer::Renderer> CreateAndInitRenderer(
    const ion::math::Range2i viewport_bounds) {
  RenderingContext context;
  context.renderer.Reset(new ion::gfx::Renderer(
      ion::gfx::GraphicsManagerPtr(new ion::gfx::GraphicsManager)));
  context.render_target_size = viewport_bounds.GetSize();
  context.shader_manager.Reset(new ion::gfxutils::ShaderManager());

  // Create renderer
  std::unique_ptr<seurat::viewer::Renderer> seurat_renderer(
      new seurat::viewer::SingleResolveRenderer());
  auto scene_factory = [viewport_bounds]() {
    return std::unique_ptr<seurat::viewer::Scene>(
        new seurat::viewer::MonoScene(viewport_bounds));
  };
  seurat_renderer->Init(scene_factory, context);
  return seurat_renderer;
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cout << "usage: " << argv[0] << " <filename>\n";
    return EXIT_FAILURE;
  }

  char* filename = argv[1];

  GLFWwindow* window = nullptr;

  if (!glfwInit()) {
    exit(EXIT_FAILURE);
  }

  glfwSetErrorCallback(ErrorCallback);

  glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_API);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);

  window = glfwCreateWindow(kWindowWidth, kWindowHeight, "Butterfly", nullptr,
                            nullptr);

  if (!window) {
    glfwTerminate();
    exit(EXIT_FAILURE);
  }

  WindowUserData user_data;
  glfwSetWindowUserPointer(window, &user_data);
  glfwSetKeyCallback(window, KeyboardCallback);
  glfwSetMouseButtonCallback(window, MouseButtonCallback);
  glfwSetCursorPosCallback(window, CursorPositionCallback);

  glfwMakeContextCurrent(window);

  std::unique_ptr<seurat::viewer::Renderer> seurat_renderer =
      CreateAndInitRenderer(
          Range2i(Point2i::Zero(), Point2i(kWindowWidth, kWindowHeight)));

  seurat_renderer->SetComponent(LoadComponent(filename));

  while (!glfwWindowShouldClose(window)) {
    seurat_renderer->RenderFrame(ComputeViewState(user_data.viewer_camera));
    glfwSwapBuffers(window);
    glfwPollEvents();
  }

  seurat_renderer.reset();

  glfwTerminate();

  return EXIT_SUCCESS;
}
