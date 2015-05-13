#ifndef _SPHERE_H_
#define _SPHERE_H_

#include <glm/glm.hpp>
#include <vector>
#include "vbo_structs.h"

// ====================================================================
// ====================================================================
// Simple implicit repesentation of a sphere, that can also be
// rasterized for use in radiosity.

class Sphere {

public:
  // CONSTRUCTOR & DESTRUCTOR
  Sphere(){
    radius = 1;
    center = glm::vec3(0,0,0);
  }

  Sphere(const glm::vec3 &c, float r) {
    center = c; radius = r;
    assert (radius >= 0); 
  }

  // used to get points along a sphere
  glm::vec3 ComputeSpherePoint(
    float s, float t, 
    const glm::vec3 center, 
    float radius);

  // for OpenGL rendering
  void setup( 
    int sphere_horiz, int sphere_vert,
    std::vector<VBOPosNormalColor> & verts, 
    std::vector<VBOIndexedTri> & tri_indices
  );

private:

  // REPRESENTATION
  glm::vec3 center;
  float radius;

};

// ====================================================================
// ====================================================================


#endif
