#version 330 core

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in vec3 vertexNormal_modelspace;
layout(location = 2) in vec3 vertexColor;

// Output data
out vec3 vertexPosition_worldspace;
out vec3 vertexNormal_worldspace;
out vec3 EyeDirection_cameraspace;
out vec3 myColor;

// Values that stay constant for the whole mesh.
uniform mat4 MVP;
uniform mat4 V;
uniform mat4 M;
uniform vec3 LightPosition_worldspace;

void main(){
  
  // Output position of the vertex, in clip space : MVP * position
  gl_Position =  MVP * vec4(vertexPosition_modelspace,1);
  
  // Position of the vertex, in worldspace : M * position
  vertexPosition_worldspace = (M * vec4(vertexPosition_modelspace,1)).xyz;
  
  // Vector that goes from the vertex to the camera, in camera space.
  // In camera space, the camera is at the origin (0,0,0).
  vec3 vertexPosition_cameraspace = ( V * M * vec4(vertexPosition_modelspace,1)).xyz;
  
  EyeDirection_cameraspace = vec3(0,0,0) - vertexPosition_cameraspace;
  
  vertexNormal_worldspace = normalize (M * vec4(vertexNormal_modelspace,0)).xyz; 
  
  // pass color to the fragment shader
  myColor = vertexColor;

}


