#version 330 core

// Interpolated values from the vertex shader
in vec3 vertexPosition_worldspace;
in vec3 EyeDirection_cameraspace;
in vec3 myColor;
in vec3 vertexNormal_worldspace;

// Ouput data
out vec3 color;

// Values that stay constant for the whole mesh.
uniform vec3 LightPosition_worldspace;


// Main method that computes color
void main(){

  vec3 LightColor = vec3(1,1,1);
  float LightPower = 4.0f; 

  // surface normal
  vec3 surface_normal =  vertexNormal_worldspace;

  // Material Properties
  vec3 MaterialDiffuseColor;

  if (myColor.x > 0.15 &&
      myColor.y > 0.15 &&
      myColor.z > 0.15 ) {

    MaterialDiffuseColor = 0.7 * vec3(1,1,1);

  } else {

    MaterialDiffuseColor = myColor;

    if( myColor.x > 0.15 &&
        myColor.y > 0.15 &&
        myColor.z > 0.15 ) {
        
        MaterialDiffuseColor *= 7.0;
      }
  }

  vec3 MaterialAmbientColor = vec3(0.3,0.3,0.3) * MaterialDiffuseColor;
  vec3 MaterialSpecularColor = vec3(0.3,0.3,0.3);

  if(!gl_FrontFacing ) {
    MaterialDiffuseColor = vec3(0.0,0.0,0.6); 
    MaterialAmbientColor = vec3(0.3,0.3,0.3) * MaterialDiffuseColor;
    MaterialSpecularColor = vec3(0.1,0.1,0.3);
    surface_normal = -surface_normal;
  }

  // Direction & distance to the light
   vec3 dirToLight = (LightPosition_worldspace - vertexPosition_worldspace);
   float distanceToLight = length( dirToLight );
   dirToLight = normalize(dirToLight);

 // float distanceToLight = length ( LightPosition_worldspace -  vertexPosition_worldspace );

  // This makes a headlamp kind of lightsource
  //vec3 dirToLight = normalize ( vec3 ( 1.0, 0.8, 2.0 ) );

  
  // Cosine of the angle between the normal and the light direction, 
  // clamped above 0
  //  - light is at the vertical of the triangle -> 1
  //  - light is perpendicular to the triangle -> 0
  //  - light is behind the triangle -> 0
  float cosTheta = clamp( dot( surface_normal,dirToLight ), 0,1 );

  // REFLECTION  
  // Eye vector (towards the camera)
  vec3 E = normalize(EyeDirection_cameraspace);
  // Direction in which the triangle reflects the light
  vec3 R = reflect(-dirToLight,surface_normal);
  // Cosine of the angle between the Eye vector and the Reflect vector,
  // clamped to 0
  //  - Looking into the reflection -> 1
  //  - Looking elsewhere -> < 1
  float cosAlpha = clamp( dot( E,R ), 0,1 );
  
  // Final assigment of color
  color = 
    MaterialAmbientColor +
    MaterialDiffuseColor * LightColor * LightPower * cosTheta / (distanceToLight) +
    MaterialSpecularColor * LightColor * LightPower * pow(cosAlpha,5) / (distanceToLight); 
  // NOTE: actually in real-world physics this should be divded by distance^2
  //    (but the dynamic range is probably too big for typical screens)
}










