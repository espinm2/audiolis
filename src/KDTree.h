#ifndef KDTREE_H
#define KDTREE_H



#include <vector>


class Particle;
class BoundingBox;

// Custom kd-tree class to store our particle system in 3d space
class KDTree{

  public:


    // Will create an balanced KDTree.
    KDTree(const std::vector<Particle *> & unsorted );

    // Used to get the children
    const Particle & rightChild(unsigned int index) const;
    const Particle &  leftChild(unsigned int index)  const;


    short unsigned int getDescriminator(unsigned int index);

    // Used for finding nearest particles
    const Particle & getNearest(
    const Particle & p, unsigned int root_index, const BoundingBox & bbox);




  private:

    // Binary Heap we will use to keep our elements
    std::vector<Particle *> elements;

};

#endif // KDTREE_H

