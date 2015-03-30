// These are nodes I will use
#ifndef BVHNODE_H
#define BVHNODE_H

#include <cassert>
#include <boundingbox.h>
#include <triangle.h>

typedef std::vector<Triangle*> TriangleVec;

// Non-Memeber Comparsion Functions comparsion functions 
bool compare_x_pos(Triangle * a, Triangle * b){
  return a->getCenter().x < b->getCenter().x; } 
bool compare_y_pos(Triangle * a, Triangle * b){
  return a->getCenter().y < b->getCenter().y; }
bool compare_z_pos(Triangle * a, Triangle * b){
  return a->getCenter().z < b->getCenter().z; }


class  BVHNode {

  public:

    // Constrcutor
    BVHNode(){}

    // Nodes will build self
    BVHNode(TriangleVec tv, int axis){

      // Sanity check
      assert(tv.size() != 0);

      if(tv.size() == 1){ // Am I a leaf

        tri_leaf = tv[0];
        left_volume = right_volume = NULL;

      } else { // Not a leaf, I am a node

        // Find bounding box of tv
        BoundingBox b;

        for(Triangle * t: tv){
          b.Extend((*t)[0]->getPos());
          b.Extend((*t)[1]->getPos());
          b.Extend((*t)[2]->getPos());
        }

        bbox = b;

        // Sort by axis
        switch(axis){

          case 0:
            std::sort(tv.begin() + i, tv.begin()+j, compare_x_pos);
            break;
          case 1:
            std::sort(tv.begin() + i, tv.begin()+j, compare_y_pos);
            break;
          case 2:
            std::sort(tv.begin() + i, tv.begin()+j, compare_z_pos);
            break;
          default:
            assert(false);
        }
  

        
        // Break up into 2 vectors
        int midpoint = tv.size() / 2;
        TriangleVec tv_left(tv.begin(), tv.begin() + midpoint);
        TriangleVec tv_right(tv.begin() + midpoint, tv.end());

        // Recurse
        left_volume  = new ( tv_left,  (axis + 1) % 3 );
        right_volume = new ( tv_right, (axis + 1) % 3 );
      
      }//else
    }//end

    // Accessors
    BVHNode * getRightVolume(){ return right_volume; }
    BVHNode * getLeftVolume(){ return left_volume; }
    BoundingBox getBoundingBox(){ return bbox; }
    bool isLeaf(){ return tri_leaf != NULL; }
    Triangle * getTriangle(){ assert(isLeaf()); return tri_leaf; }
    
  private:

    BVHNode * left_volume;    // All NULL if a leaf node
    BVHNode * right_volume;   
    BoundingBox bbox;         
    Triangle * tri_leaf; // NULL unless a leaf node
    
};

#endif
