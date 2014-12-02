#ifndef _TRIANGLE_H
#define _TRIANGLE_H

#include "boundingbox.h"
#include "edge.h"
#include <glm/glm.hpp>
#include <string>
#include "vertex.h"

// ===========================================================
// Simple half-edge data structure representation for a triangle mesh

class Triangle {

public:

  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  Triangle() {
    edge = NULL; 
    id = next_triangle_id;
    next_triangle_id++;
  }
  ~Triangle() {}

  // =========
  // ACCESSORS
  Vertex* operator[](int i) const { 
    assert (edge != NULL);
    if (i==0) return edge->getStartVertex();
    if (i==1) return edge->getNext()->getStartVertex();
    if (i==2) return edge->getNext()->getNext()->getStartVertex();
    assert(0);
  }

  Edge* getEdge() { 
    assert (edge != NULL);
    return edge; 
  }

  int getID() { return id; }
  
  std::string getMaterial() const { return mtl; }


  //get center of triangle
  glm::vec3 getCenter(){
    glm::vec3  a= edge->getStartVertex()->getPos();
    glm::vec3  b= edge->getNext()->getStartVertex()->getPos();              
    glm::vec3  c= edge->getNext()->getNext()->getStartVertex()->getPos();   
    glm::vec3 res = a + b + c;
    res = res * (1.0f/3.0f);
    return res;
}

  glm::vec3 getNormal(){
    
    glm::vec3  p1 = edge->getStartVertex()->getPos();
    glm::vec3  p2 = edge->getNext()->getStartVertex()->getPos();              
    glm::vec3  p3 = edge->getNext()->getNext()->getStartVertex()->getPos();   
  
    glm::vec3 v12 = p2;
    v12 -= p1;
    glm::vec3 v23 = p3;
    v23 -= p2;
    glm::vec3 normal = glm::normalize(glm::cross(v12,v23));
    return normal;
  
  }


  glm::vec3 getNormalRev(){

    glm::vec3  p3 = edge->getStartVertex()->getPos();
    glm::vec3  p2 = edge->getNext()->getStartVertex()->getPos();              
    glm::vec3  p1 = edge->getNext()->getNext()->getStartVertex()->getPos();   
  
    glm::vec3 v12 = p2;
    v12 -= p1;
    glm::vec3 v23 = p3;
    v23 -= p2;
    glm::vec3 normal = glm::normalize(glm::cross(v12,v23));
    return normal;
  }
  
  // =========
  // MODIFIERS
  void setEdge(Edge *e) {
    assert (edge == NULL);
    edge = e;
  }

  void setMaterial(const std::string& m){
    // Copy material
    mtl = m;
  }

protected:

  // ==============
  // REPRESENTATION
  Edge *edge;
  std::string mtl;
  int id;

  static int next_triangle_id;
};

// ===========================================================

#endif
