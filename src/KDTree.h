#ifndef KDTREE_H
#define KDTREE_H

#include <vector>
#include "particle.h"
#include "boundingbox.h"

class Particle;
class BoundingBox;
class VBOPosNormalColor;

typedef unsigned int uint;
typedef short unsigned int uint8;
typedef std::vector<Particle *>  partPtrVec;

// Custom kd-tree class to store our particle system in 3d space
class KDTree{

  public:

    // Will create an balanced KDTree.
    KDTree(){}
    void update(const partPtrVec & unsorted, const BoundingBox & bbox);
    void Optimize(uint i, uint j, partPtrVec & a, uint8 d, uint hp);

    bool ParticleSearch(const Particle * &p);

    void GatherParticles( Particle * center_particle, double gather_radius, uint heap_index, uint8 d, partPtrVec & gathered_particles, glm::vec3 minPt, glm::vec3 maxPt);

    bool Intersection(glm::vec3 min_pt, glm::vec3 max_pt, glm::vec3 sph_center, double sph_radius);
    
    // Used to get the children
    const Particle * rightChild(uint index) const;
    const Particle * leftChild(uint index)  const;

    bool hasLeft(const uint & index)const{
      return 0 <= index && index * 2 + 1 < binary_heap.size(); }

    bool hasRight(const uint & index)const{ 
      return 0 <= index && index * 2 + 2 < binary_heap.size(); }

    // Recursive function to render the tree
    void renderKDTree(uint hp, uint8 d, const glm::vec3 & minPt, const glm::vec3 & maxPt);

    // Used to rending a bounding box
    void renderBBox(const glm::vec3 &A, const glm::vec3 &B);

    void initializeVBOs();
    void setupVBOs();
    void drawVBOs();
    void cleanupVBOs();

  private:

    // Binary Heap we will use to keep our elements
    std::vector<Particle *> binary_heap;
    BoundingBox bbox;

    // Buffers for rendering
    GLuint tree_verts_VBO;
    GLuint tree_tri_indices_VBO;
    std::vector<VBOPosNormalColor> tree_verts;
    std::vector<VBOIndexedTri> tree_tri_indices;

};

#endif // KDTREE_H

