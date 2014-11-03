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
uniform int colormode;
uniform int whichshader;


//////////////////////////LIBRARY WITH NOISE//////////////////////
//
// Description : Array and textureless GLSL 2D/3D/4D simplex
// noise functions.
// Author : Ian McEwan, Ashima Arts.
// Maintainer : ijm
// Lastmod : 20110822 (ijm)
// License : Copyright (C) 2011 Ashima Arts. All rights reserved.
// Distributed under the MIT License. See LICENSE file.
// https://github.com/ashima/webgl-noise
//

vec3 mod289(vec3 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 mod289(vec4 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 permute(vec4 x) {
     return mod289(((x*34.0)+1.0)*x);
}

vec4 taylorInvSqrt(vec4 r)
{
  return 1.79284291400159 - 0.85373472095314 * r;
}

float snoise(vec3 v)
  {
  const vec2 C = vec2(1.0/6.0, 1.0/3.0) ;
  const vec4 D = vec4(0.0, 0.5, 1.0, 2.0);

// First corner
  vec3 i = floor(v + dot(v, C.yyy) );
  vec3 x0 = v - i + dot(i, C.xxx) ;

// Other corners
  vec3 g = step(x0.yzx, x0.xyz);
  vec3 l = 1.0 - g;
  vec3 i1 = min( g.xyz, l.zxy );
  vec3 i2 = max( g.xyz, l.zxy );

  // x0 = x0 - 0.0 + 0.0 * C.xxx;
  // x1 = x0 - i1 + 1.0 * C.xxx;
  // x2 = x0 - i2 + 2.0 * C.xxx;
  // x3 = x0 - 1.0 + 3.0 * C.xxx;
  vec3 x1 = x0 - i1 + C.xxx;
  vec3 x2 = x0 - i2 + C.yyy; // 2.0*C.x = 1/3 = C.y
  vec3 x3 = x0 - D.yyy; // -1.0+3.0*C.x = -0.5 = -D.y

// Permutations
  i = mod289(i);
  vec4 p = permute( permute( permute(
             i.z + vec4(0.0, i1.z, i2.z, 1.0 ))
           + i.y + vec4(0.0, i1.y, i2.y, 1.0 ))
           + i.x + vec4(0.0, i1.x, i2.x, 1.0 ));

// Gradients: 7x7 points over a square, mapped onto an octahedron.
// The ring size 17*17 = 289 is close to a multiple of 49 (49*6 = 294)
  float n_ = 0.142857142857; // 1.0/7.0
  vec3 ns = n_ * D.wyz - D.xzx;

  vec4 j = p - 49.0 * floor(p * ns.z * ns.z); // mod(p,7*7)

  vec4 x_ = floor(j * ns.z);
  vec4 y_ = floor(j - 7.0 * x_ ); // mod(j,N)

  vec4 x = x_ *ns.x + ns.yyyy;
  vec4 y = y_ *ns.x + ns.yyyy;
  vec4 h = 1.0 - abs(x) - abs(y);

  vec4 b0 = vec4( x.xy, y.xy );
  vec4 b1 = vec4( x.zw, y.zw );

  //vec4 s0 = vec4(lessThan(b0,0.0))*2.0 - 1.0;
  //vec4 s1 = vec4(lessThan(b1,0.0))*2.0 - 1.0;
  vec4 s0 = floor(b0)*2.0 + 1.0;
  vec4 s1 = floor(b1)*2.0 + 1.0;
  vec4 sh = -step(h, vec4(0.0));

  vec4 a0 = b0.xzyw + s0.xzyw*sh.xxyy ;
  vec4 a1 = b1.xzyw + s1.xzyw*sh.zzww ;

  vec3 p0 = vec3(a0.xy,h.x);
  vec3 p1 = vec3(a0.zw,h.y);
  vec3 p2 = vec3(a1.xy,h.z);
  vec3 p3 = vec3(a1.zw,h.w);

//Normalise gradients
  vec4 norm = taylorInvSqrt(vec4(dot(p0,p0), dot(p1,p1), dot(p2, p2), dot(p3,p3)));
  p0 *= norm.x;
  p1 *= norm.y;
  p2 *= norm.z;
  p3 *= norm.w;

// Mix final noise value
  vec4 m = max(0.6 - vec4(dot(x0,x0), dot(x1,x1), dot(x2,x2), dot(x3,x3)), 0.0);
  m = m * m;
  return 42.0 * dot( m*m, vec4( dot(p0,x0), dot(p1,x1),
                                dot(p2,x2), dot(p3,x3) ) );
  }
/////////////////////////END NOISE///////////////////////////////

// ----------------------------------------------
// a shader for a black & white checkerboard
vec3 checkerboard(vec3 pos) {
  // determine the parity of this point in the 3D checkerboard
  int count = 0;
  if (mod(pos.x,0.3)> 0.15) count++;
  if (mod(pos.y,0.3)> 0.15) count++;
  if (mod(pos.z,0.3)> 0.15) count++;
  if (count == 1 || count == 3) {
    return vec3(0.1,0.1,0.1);
  } else {
    return vec3(1,1,1);
  }
}


// ----------------------------------------------
// a shader for a black & white checkerboard
vec3 orange(vec3 pos, inout vec3 normal) {
  // the base color is orange!
  vec3 color = vec3(1.0,0.5,0.1);


  // FIXME: uncommenting seems to have really slow performance even
  //    when not rendering "orange".  Perhaps we should have separate
  //    shaders and not load the orange shader if it's not being used.

  // high frequency noise added to the normal for the bump map
  
  
  // ASSIGNMENT: UNCOMMENT TO SEE THE WRINKLY ORANGE!
   normal = normalize(normal+0.4*snoise(100.0*pos));

   // Tweaking the sclar value gives us large blobs
   // Smaller means bigger blobs

  
  return color;
}

// ----------------------------------------------
// a shader for a black & white checkerboard
vec3 wood(vec3 pos, inout vec3 normal) {

  //normal = normalize(normal+0.4*snoise(70.0*pos));
  // Need to convert
  float PI = 3.1415;

  // Convert to cylander
  float x = pos.x;
  float h = pos.y;
  float z = pos.z;

  // Find radius
  float r = sqrt((x*x) + (z*z));
  float angle = 0.0;

  if( x == 0.0){
    angle = PI / 2.0;
  }else{
    angle = atan(x,z);
  }

 r = r + (2 * sin( 20 *angle+h/ 150.0));
 float grain = mod(round(r), 60);

 if (grain < 40) {
    return vec3(0.6,0.4,0.2);
 }else{
    return vec3(1.0,0.5,0.1);
 }

}

vec3 watermelon(vec3 pos, inout vec3 normal) {

  normal = normalize(normal+0.4*snoise(3.0*pos));
  // Need to convert
  float PI = 3.14159265358979323846264;

  // Convert to cylander
  float x = pos.x;
  float y = pos.y;
  float z = pos.z;

  // Find radius
  float r = sqrt(  (x*x) + (y*y ) );
  float angle = 0.0;

  if( x == 0.0){
    angle = PI / 2.0;
  }else{
    angle = atan(z/y);
  }

 r = r + (2 * sin( 20 * angle +( z / 150.0)));
 float grain = mod(r, 60.0);

 if (grain < 10) {
    return vec3(0.0,0.2,0.0);
 }else{
    return vec3(0.0,0.6,0.0);
 }

}
// ----------------------------------------------
void main(){

  vec3 LightColor = vec3(1,1,1);
  float LightPower = 4.0f; 

  // surface normal
  vec3 surface_normal =  vertexNormal_worldspace;
  
  // Material properties
  vec3 MaterialDiffuseColor = myColor;
  if (whichshader == 1) {
    MaterialDiffuseColor = checkerboard(vertexPosition_worldspace);
  } else if (whichshader == 2) {
    vec3 normal2;
    MaterialDiffuseColor = wood(vertexPosition_worldspace,surface_normal);
  } else if (whichshader == 3) {
    vec3 normal3;
    MaterialDiffuseColor = watermelon(vertexPosition_worldspace,surface_normal);
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
  

  if (colormode == 0) {
    // mode 0: NO LIGHTING
    // mode 1: NO LIGHTING
    color = MaterialDiffuseColor;
  } else if (colormode == 1) {
    // mode 1: STANDARD PHONG LIGHTING (LIGHT ON)
    color = 
      MaterialAmbientColor +
      MaterialDiffuseColor * LightColor * LightPower * cosTheta / (distanceToLight) +
      MaterialSpecularColor * LightColor * LightPower * pow(cosAlpha,5) / (distanceToLight); 
    // NOTE: actually in real-world physics this should be divded by distance^2
    //    (but the dynamic range is probably too big for typical screens)
  } else if (colormode == 2) {
    // mode 2: AMBIENT ONLY (LIGHT OFF) 
    color = MaterialAmbientColor;
  }
}










