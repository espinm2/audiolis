// These are nodes I will use
#ifndef BVHNODE_H
#define BVHNODE_H

#include <cassert>
#include "boundingbox.h"
#include <algorithm>
#include "triangle.h"
#include "ray.h"
#include "hit.h"

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
    BVHNode(TriangleVec & triVec, int axis){

      TriangleVec tv = triVec; // forceful copy
      // Sanity check
      assert(tv.size() != 0);

      if(tv.size() == 1){ // Am I a leaf


        // Set my triangle
        tri_leaf = tv[0];
        left_volume = right_volume = NULL;

        BoundingBox b;
        b.Extend((*tri_leaf)[0]->getPos());
        b.Extend((*tri_leaf)[1]->getPos());
        b.Extend((*tri_leaf)[2]->getPos());
        bbox = b;

      } else { // Not a leaf, I am a node


        tri_leaf = NULL;
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
            std::sort(tv.begin(), tv.end(), compare_x_pos);
            break;
          case 1:
            std::sort(tv.begin(), tv.end(), compare_y_pos);
            break;
          case 2:
            std::sort(tv.begin(), tv.end(), compare_z_pos);
            break;
          default:
            assert(false);
        }
  

        // Break up into 2 vectors
        int midpoint = tv.size() / 2;
        TriangleVec tv_left(tv.begin(), tv.begin() + midpoint);
        TriangleVec tv_right(tv.begin() + midpoint, tv.end());

        // Recurse
        left_volume  = new BVHNode( tv_left,  (axis + 1) % 3 );
        right_volume = new BVHNode( tv_right, (axis + 1) % 3 );
      
      }//else
    }//end

    // Basic Accessors
    BVHNode * getRightVolume(){ return right_volume; }
    BVHNode * getLeftVolume(){ return left_volume; }
    BoundingBox getBoundingBox(){ return bbox; }
    bool isLeaf(){ return tri_leaf != NULL; }
    Triangle * getTriangle(){ assert(isLeaf()); return tri_leaf; }

    // Will count how many leaves I have, should match mesh triangles
    int leafCount(){

      if( isLeaf() ) { return 1; }

      return right_volume->leafCount() + left_volume->leafCount(); }

    // Collision with ray function
    bool rayHit(const Ray & r, double t1, double t2, Hit & h){
      /*std::cout << "========================================\n";
      std::cout << r  << std::endl;
      std::cout << "Interval "  << t1 << ", " << t2 << std::endl;
      std::cout << h << std::endl;*/

      if(isLeaf()){
          // Check for triangle collision
          Hit tHit;
          if(getTriangle()->rayIntersection(r,t1,t2,tHit)){
            h = tHit;
            return true;
          }else{
            return false;
          }

      }else{

        if(bbox.hitbox(r,t1,t2)){

          Hit lHit,rHit;
          bool leftHit,rightHit;

          // I know for a fact that i have a left and right
          leftHit  = left_volume->rayHit(r,t1,t2,lHit);
          rightHit = right_volume->rayHit(r,t1,t2,rHit);

          if(leftHit && rightHit){

            // Hit my left before my right
            if(lHit.getT() < rHit.getT()){

              h = lHit;

            }else{

              h = rHit;

            }

            return true;

          }else if(leftHit){

            h = lHit;
            return true;


          }else if(rightHit){

            h = rHit;
            return true;

          }else{

            return false;
          }

        }//if hit

        return false;

    }//isleaf



    }//getTriangles



    
  private:

    BVHNode * left_volume;    // All NULL if a leaf node
    BVHNode * right_volume;   
    BoundingBox bbox;         
    Triangle * tri_leaf; // NULL unless a leaf node
    
};

#endif
