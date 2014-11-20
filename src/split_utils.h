
#define EPSILON 0.0001

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <vector>
#include <cmath>
#include "ray.h"
#include "hit.h"

// This file will contain math required for particle splits

// ╔═╗╦═╗╔═╗╔╦╗╔═╗╔╦╗╦ ╦╔═╗╔═╗
// ╠═╝╠╦╝║ ║ ║ ║ ║ ║ ╚╦╝╠═╝║╣ 
// ╩  ╩╚═╚═╝ ╩ ╚═╝ ╩  ╩ ╩  ╚═╝


glm::vec3 parametric_circle_3d(glm::vec3 a, glm::vec3 b, 
    glm::vec3 c, float r, float theta);

void circle_points_on_plane( const glm::vec3 c, const glm::vec3 n, 
    const float r, const int numberPoints, std::vector<glm::vec3> &pts);

void cirlce_point_on_sphere( const glm::vec3 &center, const float radius, 
  std::vector<glm::vec3> &pts);
// Outlines a cirlce with n number of points on a plane.

// Projects these points all back onto a sphere


// ╦╔╦╗╔═╗╦  ╔═╗╔╦╗╔═╗╔╗╔╔╦╗╔═╗╔╦╗╦╔═╗╔╗╔
// ║║║║╠═╝║  ║╣ ║║║║╣ ║║║ ║ ╠═╣ ║ ║║ ║║║║
// ╩╩ ╩╩  ╩═╝╚═╝╩ ╩╚═╝╝╚╝ ╩ ╩ ╩ ╩ ╩╚═╝╝╚╝

void circle_points_on_plane( const glm::vec3 c, const glm::vec3 n, 
    const float r, const int numberPoints, std::vector<glm::vec3> &pts){

  float theta = 2 * M_PI / numberPoints;

  // Solving for a  to be orthagonal
  glm::vec3 a; a.x = 0.5; a.y = 0.5;

  if(n.z != 0){

    a.z = ( -1 * (n.x * a.x) - (n.y * a.y) )  / n.z;

  }else{

    a = glm::vec3(0,0,1);

  }

  a = glm::normalize(a);
  
  glm::vec3 b = glm::cross(n,a);
  b = glm::normalize(b);


  for(int i = 0; i < numberPoints; i++){

    // Get point
    glm::vec3 pos = parametric_circle_3d(a,b,c,r,theta*i);

    pts.push_back(pos);

  }

}

void cirlce_point_on_sphere( const glm::vec3 &center, const float radius, 
  std::vector<glm::vec3> &pts){

  // For each pts
  for( int i  = 0; i < pts.size(); i++){

    // Find direction and normalize
    glm::vec3 dir = glm::normalize(pts[i]- center);
    
    // Find new location of point
    pts[i] = (dir * radius)  + center;
  }
}

glm::vec3 parametric_circle_3d(glm::vec3 a, glm::vec3 b, 
    glm::vec3 c, float r, float theta){

  // Make sure angle given is within range and you have directions
  assert( 0 <= theta &&  theta <= 2*M_PI);
  // assert( fabs(glm::length(a) - 1) < EPSILON );
  // assert( fabs(glm::length(b) - 1) < EPSILON );

  float x = c.x + r * cos(theta) * a.x + r * sin(theta) * b.x;
  float y = c.y + r * cos(theta) * a.y + r * sin(theta) * b.y;
  float z = c.z + r * cos(theta) * a.z + r * sin(theta) * b.z;

  return glm::vec3(x,y,z);

}
