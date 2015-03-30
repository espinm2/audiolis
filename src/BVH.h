#ifndef BVH_H
#define BVH_H

typedef std::vector<Triangle*> triangleVec;

class Mesh;
class BVHNode;
class Triangle;


class BVH{

  public:

    // Constructor
    BVH(Mesh * mesh);
    void BVHaux( BVHNode * node, TriangleVec, int axis);

  private:

    // Root of our tree
    BVHNode * root;

}

#endif
