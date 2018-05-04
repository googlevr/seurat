# Copyright 2017 Google Inc. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS-IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Generates a JSON manifest and a Maya camera rig for Seurat.

   Example usage:
     CreateRig(headbox_min=[-0.5, -0.5, -0.5],
               headbox_max=[0.5, 0.5, 0.5],
               num_view_groups=16,  # Should be a power of two.
               image_size=1024,
               near_clip=0.1,
               far_clip=100.0,
               depth_type='EYE_Z',
               depth_channel_name='A',
               color_file_path_pattern='%s_color.%04d.exr',
               depth_file_path_pattern='%s_depth.%04d.exr',
               json_file_path='./manifest.json')
"""
import json
import math
import operator


def ProjectPoint(matrix, point):
  """Projects a 3D point using a 4x4 matrix.

  Args:
    matrix: A 4x4 matrix represented as a list of 16 floats.
    point: A 3D point represented as a list of 3 floats.

  Returns:
    The projected point, represented as a list of 3 floats.
  """
  result_hom = [0.0, 0.0, 0.0, 0.0]
  for row in xrange(4):
    for col in xrange(3):
      result_hom[row] += matrix[4 * row + col] * point[col]
    # point.w = 1.0 implicitly
    result_hom[row] += matrix[4 * row + 3]
  w = result_hom[3]
  return map(operator.div, result_hom[0:3], [w, w, w])


def WorldFromEyeMatrixFromFace(face_name):
  """Creates world-from-eye matrix for the given face of a cube map.

  Args:
    face_name: Name of the face. Must be one of 'front', 'back', 'left',
        'right', 'bottom', 'top'.

  Returns:
    The world-from-eye matrix for the given face as a list in row-major order.

  Raises:
    ValueError: face_name is not the name of a cube map face.
  """
  # pylint: disable=bad-whitespace
  # pylint: disable=bad-continuation
  if face_name is 'front':
    return [ 1.0,  0.0,  0.0,  0.0,
             0.0,  1.0,  0.0,  0.0,
             0.0,  0.0,  1.0,  0.0,
             0.0,  0.0,  0.0,  1.0]  # pyformat: disable
  elif face_name is 'back':
    return [-1.0,  0.0,  0.0,  0.0,
             0.0,  1.0,  0.0,  0.0,
             0.0,  0.0, -1.0,  0.0,
             0.0,  0.0,  0.0,  1.0]  # pyformat: disable
  elif face_name is 'left':
    return [ 0.0,  0.0,  1.0,  0.0,
             0.0,  1.0,  0.0,  0.0,
            -1.0,  0.0,  0.0,  0.0,
             0.0,  0.0,  0.0,  1.0]  # pyformat: disable
  elif face_name is 'right':
    return [ 0.0,  0.0, -1.0,  0.0,
             0.0,  1.0,  0.0,  0.0,
             1.0,  0.0,  0.0,  0.0,
             0.0,  0.0,  0.0,  1.0]  # pyformat: disable
  elif face_name is 'bottom':
    return [ 1.0,  0.0,  0.0,  0.0,
             0.0,  0.0,  1.0,  0.0,
             0.0, -1.0,  0.0,  0.0,
             0.0,  0.0,  0.0,  1.0]  # pyformat: disable
  elif face_name is 'top':
    return [ 1.0,  0.0,  0.0,  0.0,
             0.0,  0.0, -1.0,  0.0,
             0.0,  1.0,  0.0,  0.0,
             0.0,  0.0,  0.0,  1.0]  # pyformat: disable
  else:
    raise ValueError('Invalid face_name')


def CubeFaceProjectionMatrix(near, far):
  """Creates a cube-face 90 degree FOV projection matrix.

  The created matrix is an OpenGL-style projection matrix.

  Args:
    near: Eye-space Z position of the near clipping plane.
    far: Eye-space Z position of the far clipping plane.

  Returns:
    The clip-from-eye matrix as a list in row-major order.

  Raises:
    ValueError: Invalid clip planes. near <= 0.0 or far <= near.
  """
  if near <= 0.0:
    raise ValueError('near must be positive.')

  if far <= near:
    raise ValueError('far must be greater than near.')

  left = -near
  right = near
  bottom = -near
  top = near
  a = (2.0 * near) / (right - left)
  b = (2.0 * near) / (top - bottom)
  c = (right + left) / (right - left)
  d = (top + bottom) / (top - bottom)
  e = (near + far) / (near - far)
  f = (2.0 * near * far) / (near - far)

  # pylint: disable=bad-whitespace
  return [a,   0.0,  c,   0.0,
          0.0, b,    d,   0.0,
          0.0, 0.0,  e,   f,
          0.0, 0.0, -1.0, 0.0]  # pyformat: disable


def RadicalInverse(a, base):
  """Computes the radical inverse of |a| in base |base|.

  Args:
    a: The integer number for which the radical inverse is computed.
    base: The radical inverse is computed in this base (integer).

  Returns:
    The radical inverse as a float in the range [0.0, 1.0).
  """
  reversed_digits = 0
  base_n = 1
  # Compute the reversed digits, base b.
  while a > 0:
    next_a = a / base
    digit = a - next_a * base
    reversed_digits = reversed_digits * base + digit
    base_n *= base
    a = next_a
  # Only when done are the reversed digits divided by b^n.
  return min(reversed_digits / float(base_n), 1.0)


def PointInBox(box_min, box_max, sample):
  """Computes a sample point inside a box with arbitrary number of dimensions.

  Args:
    box_min: A list of floats representing the lower bounds of the box.
    box_max: A list of floats representing the upper bounds of the box.
    sample: A list of floats in the range [0.0, 1.0] representing the
      relative sample position in the box.

  Returns:
    A list of floats, representing the absolute position of the sample in
    the box.
  """
  delta = map(operator.sub, box_max, box_min)
  offset = map(operator.mul, delta, sample)
  position = map(operator.add, box_min, offset)
  return position


def Distance(point_a, point_b):
  """Computes the euclidean distance between two points.

  The points can have an aribtrary number of dimensions.

  Args:
    point_a: A list of numbers representing the first point.
    point_b: A list of numbers representing the second point.

  Returns:
    The euclidean distance as a float.
  """
  delta = map(operator.sub, point_a, point_b)
  delta_sqr = map(operator.mul, delta, delta)
  distance_sqr = 0.0
  for element in delta_sqr:
    distance_sqr += element
  return math.sqrt(distance_sqr)


def RotateCamera(camera_name, face_name):
  """Rotates a Maya camera node to look at a given cube map face.

  Args:
    camera_name: Name of the Maya camera's transform node.
    face_name: Name of the cube map face.

  Raises:
    ValueError: face is not a valid cube map face name.
  """
  # Disable the undefined-variable lint error, because the Maya package is not
  # defined in the environment where the linter runs.
  #
  # pylint: disable=undefined-variable
  if face_name is 'front':
    pass
  elif face_name is 'back':
    maya.cmds.setAttr(camera_name + '.rotateY', 180)
  elif face_name is 'left':
    maya.cmds.setAttr(camera_name + '.rotateY', 90)
  elif face_name is 'right':
    maya.cmds.setAttr(camera_name + '.rotateY', -90)
  elif face_name is 'bottom':
    maya.cmds.setAttr(camera_name + '.rotateX', -90)
  elif face_name is 'top':
    maya.cmds.setAttr(camera_name + '.rotateX', 90)
  else:
    raise ValueError('Invalid face_name')


def GenerateCameraPositions(headbox_min, headbox_max, num_cameras):
  """Generates camera positions in a headbox.

  Camera posittions are computed as a 3D Hammersley point set. The points are
  transformed such that their bounding box is exactly equal to the headbox. The
  points are then sorted according to distance to the headbox center. Finally,
  the point that is closest to the headbox center is replaced by the headbox
  center itself to include a view from the reference camera.

  Args:
    headbox_min: The lower bounds of the headbox as a list of 3 floats.
    headbox_max: The upper bounds of the headbox as a list of 3 floats.
    num_cameras: The number of cameras to generate. Should be a power of two.

  Returns:
    A list of 3D points (each a list of 3 floats), representing the positions
    of the generated cameras.

  Raises:
    ValueError: num_cameras is not positive.

  """
  if num_cameras <= 0:
    raise ValueError('num_cameras must be positive')

  if num_cameras == 1:
    # Use the headbox center if a single camera position is requested.
    return [PointInBox(headbox_min, headbox_max, [0.5, 0.5, 0.5])]

  samples = []
  max_sample = [0.0, 0.0, 0.0]
  for i in xrange(num_cameras):
    # Use a 3D Hammersley point set for the samples.
    sample = [
        i / float(num_cameras),
        RadicalInverse(i, 2),
        RadicalInverse(i, 3)
    ]
    for dim in xrange(3):
      max_sample[dim] = max(max_sample[dim], sample[dim])
    samples.append(sample)

  headbox_center = PointInBox(headbox_min, headbox_max, [0.5, 0.5, 0.5])
  camera_positions = []

  for sample in samples:
    # Normalize the samples so that their bounding box is the unit cube.
    for dim in xrange(3):
      sample[dim] /= max_sample[dim]
    position = PointInBox(headbox_min, headbox_max, sample)
    camera_positions.append(position)

  sorted_positions = sorted(
      camera_positions, key=lambda point: Distance(point, headbox_center))
  # Replace the point closest to the headbox center by the headbox center
  # itself.
  sorted_positions[0] = PointInBox(headbox_min, headbox_max, [0.5, 0.5, 0.5])
  return sorted_positions


def CreateCameras(camera_positions, near_clip, far_clip):
  """Creates and animates the Maya cameras for the rig.

  Six cameras, one for each cube face, are generated. Each camera is configured
  with a square viewport and the given near and far clipping planes. This method
  also adjusts the Maya timeline to exactly contain the frames for the rig
  animation. Each of the six cameras will get one keyframe per camera position.

  Args:
    camera_positions: A list of 3D points (each a list of 3 floats) representing
        the positions of the cameras.
    near_clip: Eye-space Z position of the near clipping planes.
    far_clip: Eye-space Z position of the far clipping planes.
  """
  # Disable the undefined-variable lint error, because the Maya package is not
  # defined in the environment where the linter runs.
  #
  # pylint: disable=undefined-variable
  start_time = 0
  end_time = len(camera_positions) - 1
  maya.cmds.playbackOptions(
      animationStartTime=start_time,
      animationEndTime=end_time,
      minTime=start_time,
      maxTime=end_time)
  for face in ['front', 'back', 'left', 'right', 'bottom', 'top']:
    # Create a cube face camera and rotate it.
    camera_name = maya.cmds.camera(
        name='seurat_' + face,
        focalLength=12.7,
        horizontalFilmAperture=1,
        verticalFilmAperture=1,
        nearClipPlane=near_clip,
        farClipPlane=far_clip)[0]
    RotateCamera(camera_name, face)

    # Set translation keyframes for all positions on this camera.
    for view_group_index, position in enumerate(camera_positions):
      maya.cmds.setKeyframe(
          camera_name, at='translateX', t=view_group_index, v=position[0])
      maya.cmds.setKeyframe(
          camera_name, at='translateY', t=view_group_index, v=position[1])
      maya.cmds.setKeyframe(
          camera_name, at='translateZ', t=view_group_index, v=position[2])


def CreateViewGroups(headbox_center, camera_positions, image_size, near_clip,
                     far_clip, depth_type, depth_channel_name,
                     color_file_path_pattern, depth_file_path_pattern):
  """Creates and returns the view groups for the JSON output.

  Args:
    headbox_center: Center of the headbox as a list of 3 floats.
    camera_positions: Positions of the cameras as a list of 3D points (each a
        list of 3 floats).
    image_size: Size of the output images in pixels.
    near_clip: Eye-space Z position of the near clipping planes.
    far_clip: Eye-space Z position of the far clipping planes.
    depth_type: A string representing the depth encoding. Valid values are:
        'WINDOW_Z' (window-space Z coordinate in the range [0.0, 1.0]),
        'EYE_Z' (negated eye-space Z coordinate in the range [0.0, inf);
                 Arnold's encoding),
        'RAY_DEPTH' (distance to eye).
    depth_channel_name: Name of the depth channel in the output file.
        Commonly used values are 'R' (VRay) and 'A' (Arnold).
    color_file_path_pattern: File name pattern for color images. Must
        contain a placeholder for a string (face name) and an integer (view
        group number). Example: '%s.%04d.exr' for file names
        'front.0000.exr', 'front.0001.exr', ... , 'top.9999.exr'.
    depth_file_path_pattern: File name pattern for depth images. Must
        contain a placeholder for a string (face name) and an integer (view
        group number). Example: '%s.%04d.exr' for file names
        'front.0000.exr', 'front.0001.exr', ... , 'top.9999.exr'.

  Returns:
    A dictionary representing the view groups.
  """
  view_groups = []
  for view_group_index, absolute_position in enumerate(camera_positions):
    views = []
    for face in ['front', 'back', 'left', 'right', 'bottom', 'top']:
      # Camera position relative to headbox center.
      position = map(operator.sub, absolute_position, headbox_center)
      clip_from_eye_matrix = CubeFaceProjectionMatrix(near_clip, far_clip)
      world_from_eye_matrix = WorldFromEyeMatrixFromFace(face)
      # Set translation component of world-from-eye matrix.
      for i in xrange(3):
        world_from_eye_matrix[4 * i + 3] = position[i]

      # Create camera object
      camera = {
          'image_width': image_size,
          'image_height': image_size,
          'clip_from_eye_matrix': clip_from_eye_matrix,
          'world_from_eye_matrix': world_from_eye_matrix,
          'depth_type': depth_type
      }

      # Create view object and add it to the view groups
      color_image_path = (color_file_path_pattern % (face, view_group_index))
      depth_image_path = (depth_file_path_pattern % (face, view_group_index))
      view = {
          'projective_camera': camera,
          'depth_image_file': {
              'color': {
                  'path': color_image_path,
                  'channel_0': 'R',
                  'channel_1': 'G',
                  'channel_2': 'B',
                  'channel_alpha': 'A'
              },
              'depth': {
                  'path': depth_image_path,
                  'channel_0': depth_channel_name
              }
          }
      }
      views.append(view)
    view_group = {'views': views}
    view_groups.append(view_group)

  # Return the view_groups as a Python list.
  return view_groups


def CreateRig(headbox_min,
              headbox_max,
              num_view_groups,
              image_size,
              near_clip,
              far_clip,
              depth_type,
              depth_channel_name,
              color_file_path_pattern,
              depth_file_path_pattern,
              json_file_path,
              json_only=False):
  """Creates a Maya camera rig and JSON manifest for Seurat.

  Args:
    headbox_min: List of three floats representing the lower bounds of the
        headbox in world-space.
    headbox_max: List of three floats representing the upper bounds of the
        headbox in world-space.
    num_view_groups: Number of view groups (camera positions) to generate.
      Must be a power of two.
    image_size: Resolution of the output images in pixels.
    near_clip: Eye-space Z position of the near clipping planes.
    far_clip: Eye-space Z position of the far clipping planes.
    depth_type: A string representing the depth encoding. Valid values are:
        'WINDOW_Z' (window-space Z coordinate in the range [0.0, 1.0]),
        'EYE_Z' (negated eye-space Z coordinate in the range [0.0, inf);
                 Arnold's encoding),
        'RAY_DEPTH' (distance to eye).
    depth_channel_name: Name of the depth channel in the output file.
        Commonly used values are 'R' (VRay) and 'A' (Arnold).
    color_file_path_pattern: File name pattern for color images. Must
        contain a placeholder for a string (face name) and an integer (view
        group number). Example: '%s.%04d.exr' for file names
        'front.0000.exr', 'front.0001.exr', ... , 'top.9999.exr'.
    depth_file_path_pattern: File name pattern for depth images. Must
        contain a placeholder for a string (face name) and an integer (view
        group number). Example: '%s.%04d.exr' for file names
        'front.0000.exr', 'front.0001.exr', ... , 'top.9999.exr'.
    json_file_path: Path to the output JSON manifest file.
    json_only: A boolean value. If true, the Maya camera generation step is
        bypassed.
  """
  # Compute the positions of the cameras.
  camera_positions = GenerateCameraPositions(headbox_min, headbox_max,
                                             num_view_groups)

  # Generate the six Maya cameras and keyframe their positions.
  if not json_only:
    CreateCameras(camera_positions, near_clip, far_clip)

  # Compute the headbox center.
  headbox_center = PointInBox(headbox_min, headbox_max, [0.5, 0.5, 0.5])

  # Generate the JSON manifest and write it to the file.
  view_groups = CreateViewGroups(headbox_center, camera_positions, image_size,
                                 near_clip, far_clip, depth_type,
                                 depth_channel_name, color_file_path_pattern,
                                 depth_file_path_pattern)
  json_string = json.dumps({'view_groups': view_groups}, indent=2)
  with open(json_file_path, 'w') as json_file:
    json_file.write(json_string)
