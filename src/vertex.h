#ifndef _VERTEX_H
#define _VERTEX_H

#include <cstdlib>
#include <glm/glm.hpp>

// ==========================================================

class Vertex {

public:

  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  Vertex(int i, const glm::vec3 &pos) : position(pos), index(i) {}
  ~Vertex() {}
  
  // =========
  // ACCESSORS
  int getIndex() const { return index; }
  double x() const { return position.x; }
  double y() const { return position.y; }
  double z() const { return position.z; }
  const glm::vec3& getPos() const { return position; }
  const glm::vec3& getGouraudNormal() const { return gouraud_normal; }

  // =========
  // MODIFIERS
  void setPos(const glm::vec3 &v) { position = v; }
  void clearGouraudNormal() { gouraud_normal = glm::vec3(0,0,0); }
  void incrGouraudNormal(const glm::vec3 &v) { gouraud_normal += v; }
  void normalizeGouraudNormal() { 
    gouraud_normal = glm::normalize(gouraud_normal); }

private:

  // ==============
  // REPRESENTATION
  glm::vec3 position;
  glm::vec3 gouraud_normal;

  // this is the index from the original .obj file.
  // technically not part of the half-edge data structure, 
  // but we use it for hashing
  int index;  

  // NOTE: the vertices don't know anything about adjacency.  In some
  // versions of this data structure they have a pointer to one of
  // their incoming edges.  However, this data is very complicated to
  // maintain during mesh manipulation, so it has been omitted.

};

// ==========================================================

#endif

