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
    void update(const std::vector<Particle *> & unsorted, const BoundingBox & bbox);
    void Optimize(uint i, uint j, partPtrVec & a, uint8 d, uint hp);
    
    // Used to get the children
    const Particle * rightChild(uint index) const;
    const Particle * leftChild(uint index)  const;

    // Used to set children
    void rightChild(uint index, Particle * p);
    void leftChild(uint index, Particle * p);

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

