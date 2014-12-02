#include "geometry_utils.h"
#include "argparser.h" // This if only for the MTRand obj
#include <cstdlib>

// ╦╔╦╗╔═╗╦  ╔═╗╔╦╗╔═╗╔╗╔╔╦╗╔═╗╔╦╗╦╔═╗╔╗╔
// ║║║║╠═╝║  ║╣ ║║║║╣ ║║║ ║ ╠═╣ ║ ║║ ║║║║
// ╩╩ ╩╩  ╩═╝╚═╝╩ ╩╚═╝╝╚╝ ╩ ╩ ╩ ╩ ╩╚═╝╝╚╝

inline float DistanceBetweenTwoPoints(const glm::vec3 &p1, 
    const glm::vec3 &p2) {

  glm::vec3 v = p1-p2;
  return glm::length(v);
}

inline float AreaOfTriangle(float a, float b, float c) {
  // from the lengths of the 3 edges, compute the area
  // Area of Triangle = (using Heron's Formula)
  //  sqrt[s*(s-a)*(s-b)*(s-c)]
  //    where s = (a+b+c)/2
  // also... Area of Triangle = 0.5 * x * c
  float s = (a+b+c) / (float)2;
  return sqrt(s*(s-a)*(s-b)*(s-c));
}

inline float AreaOfTriangle(const glm::vec3 &a, 
    const glm::vec3 &b, const glm::vec3 &c) {

  float aside = DistanceBetweenTwoPoints(a,b);
  float bside = DistanceBetweenTwoPoints(b,c);
  float cside = DistanceBetweenTwoPoints(c,a);
  return AreaOfTriangle(aside,bside,cside);
}

// compute the perfect mirror direction
glm::vec3 MirrorDirection(const glm::vec3 &normal, 
    const glm::vec3 &incoming) {

  float dot = glm::dot(incoming,normal);
  glm::vec3 r = (incoming*-1.0f) + normal * (2 * dot);
  return r*-1.0f;

}

void circle_points_on_plane( const glm::vec3 c, const glm::vec3 n, 
    const float r, const int numberPoints, 
    std::vector<glm::vec3> &pts, ArgParser * &args){

  float theta = 2 * M_PI / numberPoints;

  // Solving for a  to be orthagonal // should be random
  glm::vec3 a; a.x = 0.5; a.y = 0.5;

  //glm::vec3 a; 
  //a.x =  args->randomGen.rand();
  //a.y =  args->randomGen.rand();
  // std::cout << "x: " << a.x << " y: " << a.y << std::endl;

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
  for(unsigned int i  = 0; i < pts.size(); i++){

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

inline glm::vec3 compute_normal(const glm::vec3 &p1, const glm::vec3 &p2, 
    const glm::vec3 &p3){
  glm::vec3 v12 = p2;
  v12 -= p1;
  glm::vec3 v23 = p3;
  v23 -= p2;
  glm::vec3 normal = glm::normalize(glm::cross(v12,v23));
  return normal;
}


bool plane_intersect(
    const Ray &r, Hit &h, 
    const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c, 
    bool intersect_backfacing){

  // insert the explicit equation for the 
  // ray into the implicit equation of the plane

  // equation for a plane
  // ax + by + cz = d;
  // normal . p + direction = 0
  // plug in ray
  // origin + direction * t = p(t)
  // origin . normal + t * direction . normal = d;
  // t = d - origin.normal / direction.normal;

  // Get normal
  glm::vec3 normal = compute_normal(a,b,c);

  float d = glm::dot(normal,a);

  float numer = d - glm::dot(r.getOrigin(),normal);
  float denom = glm::dot(r.getDirection(),normal);

  if (denom == 0) return 0;  // parallel to plane

  if (!intersect_backfacing && glm::dot(normal,r.getDirection()) >= 0) 
    return 0; // hit the backside

  double t = numer / denom;

  if (t > EPSILON && t < h.getT()) {
    h.set(t,normal);
    assert (h.getT() >= EPSILON);
    return 1;
  }

  return 0;
}


bool triangle_intersect( 
    const Ray &r, Hit &h, 
    const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c, 
    bool intersect_backfacing){

  // compute the intersection with the plane of the triangle
  Hit h2 = Hit(h);

  if (!plane_intersect(r,h2,a,b,c,intersect_backfacing)) return 0;  

  // figure out the barycentric coordinates:
  glm::vec3 Ro = r.getOrigin();
  glm::vec3 Rd = r.getDirection();
  // [ ax-bx   ax-cx  Rdx ][ beta  ]     [ ax-Rox ] 
  // [ ay-by   ay-cy  Rdy ][ gamma ]  =  [ ay-Roy ] 
  // [ az-bz   az-cz  Rdz ][ t     ]     [ az-Roz ] 
  
  // solve for beta, gamma, & t using Cramer's rule
  glm::mat3 detA_mat(a.x-b.x, a.x-c.x, Rd.x,
                     a.y-b.y, a.y-c.y, Rd.y,
                     a.z-b.z, a.z-c.z, Rd.z);

  float detA = glm::determinant(detA_mat);

  if (fabs(detA) <= 0.000001) return 0;

  assert (fabs(detA) >= 0.000001);
  
  glm::mat3 beta_mat(a.x-Ro.x, a.x-c.x, Rd.x,
                     a.y-Ro.y, a.y-c.y, Rd.y,
                     a.z-Ro.z, a.z-c.z, Rd.z);
  
  glm::mat3 gamma_mat(a.x-b.x, a.x-Ro.x, Rd.x,
                      a.y-b.y, a.y-Ro.y, Rd.y,
                      a.z-b.z, a.z-Ro.z, Rd.z);


  float beta = glm::determinant(beta_mat) / detA;
  float gamma = glm::determinant(gamma_mat) / detA;

  if (beta >= -0.00001 && beta <= 1.00001 &&
      gamma >= -0.00001 && gamma <= 1.00001 &&
      beta + gamma <= 1.00001) {
    h = h2;
    // interpolate the texture coordinates //
    // float alpha = 1 - beta - gamma;
    // float t_s = alpha * a->get_s() + beta * b->get_s() + gamma * c->get_s();
    // float t_t = alpha * a->get_t() + beta * b->get_t() + gamma * c->get_t();
    // h.setTextureCoords(t_s,t_t);
    assert (h.getT() >= EPSILON);
    return 1;
  }

  return 0;
}
