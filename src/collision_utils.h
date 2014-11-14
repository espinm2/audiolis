#define EPSILON 0.0001

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include "ray.h"
#include "hit.h"

// This file will contain math required for collision detection
// ╔═╗╦═╗╔═╗╔╦╗╔═╗╔╦╗╦ ╦╔═╗╔═╗
// ╠═╝╠╦╝║ ║ ║ ║ ║ ║ ╚╦╝╠═╝║╣ 
// ╩  ╩╚═╚═╝ ╩ ╚═╝ ╩  ╩ ╩  ╚═╝

inline glm::vec3 compute_normal(const glm::vec3 &p1, const glm::vec3 &p2, 
    const glm::vec3 &p3);

bool plane_intersect(
    const Ray &r, Hit &h, 
    const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c, 
    bool intersect_backfacing);

bool triangle_intersect( 
    const Ray &r, Hit &h, 
    const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c, 
    bool intersect_backfacing);

// ╦╔╦╗╔═╗╦  ╔═╗╔╦╗╔═╗╔╗╔╔╦╗╔═╗╔╦╗╦╔═╗╔╗╔
// ║║║║╠═╝║  ║╣ ║║║║╣ ║║║ ║ ╠═╣ ║ ║║ ║║║║
// ╩╩ ╩╩  ╩═╝╚═╝╩ ╩╚═╝╝╚╝ ╩ ╩ ╩ ╩ ╩╚═╝╝╚╝

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

  float t = numer / denom;

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
