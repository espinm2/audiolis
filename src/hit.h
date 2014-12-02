#ifndef _HIT_H_
#define _HIT_H_

#include <glm/glm.hpp>
#include <float.h>
#include <ostream>
#include <string>

#include "ray.h"

// Hit class mostly copied from Peter Shirley and Keith Morley
// ====================================================================
// ====================================================================

class Hit {
  
public:

  // CONSTRUCTOR & DESTRUCTOR
  Hit() { 
    t = FLT_MAX;
    normal = glm::vec3(0,0,0); 
  }
  Hit(const Hit &h) { 
    t = h.t; 
    normal = h.normal; 
    mtlHit = h.mtlHit;
    
  }
  ~Hit() {}

  // ACCESSORS
  float getT() const { return t; }
  glm::vec3 getNormal() const { return normal; }
  std::string getMaterial() const { return mtlHit; }

  // MODIFIER
  void set(float _t, glm::vec3 n) {
    t = _t; normal = n; 
  }

  void setMaterial(const std::string & mat){
    mtlHit = mat;
  }

private: 

  // REPRESENTATION
  float t;
  glm::vec3 normal;
  std::string mtlHit;

};

inline std::ostream &operator<<(std::ostream &os, const Hit &h) {
  os << "Hit <" << h.getT() << ", < "
    << h.getNormal().x << "," 
    << h.getNormal().y << "," 
    << h.getNormal().z << " > > ";
  return os;
}
// ====================================================================
// ====================================================================

#endif
