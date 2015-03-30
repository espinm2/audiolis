// Implementation of bounding volume heirachy
#include "mesh.h"   
#include "triangle.h"
#include "BVHNode.h"
#include "BVH.h"

typedef std::vector<Triangle*> TriangleVec;

// Constructor that calls aux
BVH::BVH(Mesh * mesh){

  // Push back all our triangles into a vector for easy sorting
  TriangleVec tv;
  for ( triangleshashtype::iterator iter = mesh->triangles.begin();
        iter != mesh->triangles.end(); iter++) {
      tv.push_back(iter->second);
  }


  BoundingBox bbox = mesh->getBoundingBox();

  root = new BVHNode(bbox,)



}


// Recursivly builds the tree
void BVH::BVHaux( BVHNode * node, TriangleVec tv, int axis){
  // Create left and right
}

