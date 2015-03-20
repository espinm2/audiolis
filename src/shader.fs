#version 330 core

// Interpolated values from the vertex shader
in vec3 vertexPosition_worldspace;
in vec3 EyeDirection_cameraspace;
in vec3 vertexNormal_worldspace;
in vec4 myColor;

// Ouput data
out vec4 color;

// Values that stay constant for the whole mesh.
uniform vec3 LightPosition_worldspace;
uniform int whichshader;

// ----------------------------------------------
void main(){

  vec3 LightColor = vec3(1,1,1);
  float LightPower = 4.0f; 

  // surface normal
  vec3 surface_normal =  vertexNormal_worldspace;
  
  // Get materials defuse color
  
  // Material Diffuse
  vec3 MaterialDiffuseColor = vec3(0.5,0.5,0.5);


  // SHADER FOR TRIANGULATED MESH /////////////////////////////////////////////
  if ( whichshader == 0)
  { 

    /*
    // REMOVE ME IF YOU DO NOT WANT WIREFRAME
    // Material properties
      vec3 MaterialDiffuseColor;
      if ( (wireframe == 0) ||
           ( myColor.x > 0.1 &&
             myColor.y > 0.1 &&
             myColor.z > 0.1) ) {
        MaterialDiffuseColor = 0.7 * vec3(1,1,1);
      } else {
        MaterialDiffuseColor = mycolor;
        if (mycolor.x < 0.1 &&
            mycolor.y < 0.1 &&
            mycolor.z < 0.1) {
          MaterialDiffuseColor *= 7.0; 
        }
      }
      // End wireframe
      */

      MaterialDiffuseColor.r = myColor.r;
      MaterialDiffuseColor.g = myColor.g;
      MaterialDiffuseColor.b = myColor.b;

      vec3 MaterialAmbientColor = vec3(0.3,0.3,0.3) * MaterialDiffuseColor;
      vec3 MaterialSpecularColor = vec3(0.3,0.3,0.3);
      if(!gl_FrontFacing ) 
      {
        MaterialDiffuseColor = vec3(0.0,0.0,0.6); 
        MaterialAmbientColor = vec3(0.3,0.3,0.3) * MaterialDiffuseColor;
        MaterialSpecularColor = vec3(0.1,0.1,0.3);
        surface_normal = -surface_normal;
      }

      // Direction & distance to the light
      vec3 dirToLight = (LightPosition_worldspace - vertexPosition_worldspace);
      float distanceToLight = length( dirToLight );
      dirToLight = normalize(dirToLight);
      

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
      

      vec3 no_alpha_color = 
        MaterialAmbientColor +
        MaterialDiffuseColor * LightColor * LightPower * cosTheta / (distanceToLight) +
        MaterialSpecularColor * LightColor * LightPower * pow(cosAlpha,5) / (distanceToLight); 

        color.r = no_alpha_color.r;
        color.g = no_alpha_color.g;
        color.b = no_alpha_color.b;
        color.a = myColor.a;
  }

  // SHADER FOR GLPOINTS //////////////////////////////////////////////////////
  else if (whichshader == 1) 
  {
      // Flat colors
      color = myColor;
  }

}//main










