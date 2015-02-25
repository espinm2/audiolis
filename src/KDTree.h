#ifndef KDTREE_H
#define KDTREE_H

#include <vector>

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
    KDTree(const std::vector<Particle *> & unsorted );
    void Optimize(uint i, uint j, partPtrVec & a, uint8 d, uint hp);
    
    // Used to get the children
    const Particle * rightChild(uint index) const;
    const Particle * leftChild(uint index)  const;

    // Used to set children
    void rightChild(uint index, Particle * p);
    void leftChild(uint index, Particle * p);

    // Used to render  for debuging
    void renderKDTree(std::vector<VBOPosNormalColor> & outline_verts) const;

  private:

    // Binary Heap we will use to keep our elements
    std::vector<Particle *> binary_heap;

};

#endif // KDTREE_H

