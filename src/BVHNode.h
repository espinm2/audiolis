// These are nodes I will use
#ifndef BVHNODE_H
#define BVHNODE_H

#include <cassert>
#include "boundingbox.h"
#include <algorithm>
#include "triangle.h"
#include "ray.h"

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

    // Accessors
    BVHNode * getRightVolume(){ return right_volume; }
    BVHNode * getLeftVolume(){ return left_volume; }
    BoundingBox getBoundingBox(){ return bbox; }
    bool isLeaf(){ return tri_leaf != NULL; }
    Triangle * getTriangle(){ assert(isLeaf()); return tri_leaf; }

    int leafCount(){

      if( isLeaf() ) { return 1; }

      return right_volume->leafCount() + left_volume->leafCount();
      
    }

    // Collision with ray function
    void getTriangles(Ray & r,  double time_step, TriangleVec & tv ){
      // If I am a leaf take me into being concidered for intersection
      // if( isLeaf() ) { tv.push_back(getTriangle()); return; } 

      // Check to see if this ray hits me
      glm::vec3 d = r.getDirection();
      glm::vec3 e = r.getOrigin();

      // Intersection calculations
      double tx_min, tx_max, // Assign of these
             ty_min, ty_max,
             tz_min, tz_max,
             a;

      // X's
      a = 1.0 / d.x;
      if(a >= 0 ){

        tx_min =  a * ( bbox.getMin().x - e.x );
        tx_max =  a * ( bbox.getMax().x - e.x );

      }else{

        tx_min =  a * ( bbox.getMax().x - e.x );
        tx_max =  a * ( bbox.getMin().x - e.x );
      
      }

      // Y's
      a = 1.0 / d.y;
      if(a >= 0 ){

        ty_min =  a * ( bbox.getMin().y - e.y );
        ty_max =  a * ( bbox.getMax().y - e.y );

      }else{

        ty_min =  a * ( bbox.getMax().y - e.y );
        ty_max =  a * ( bbox.getMin().y - e.y );
      
      }


      // Z's
      a = 1.0 / d.z;
      if(a >= 0 ){

        tz_min =  a * ( bbox.getMin().z - e.z );
        tz_max =  a * ( bbox.getMax().z - e.z );

      }else{

        tz_min =  a * ( bbox.getMax().z - e.z );
        tz_max =  a * ( bbox.getMin().z - e.z );
      
      }

      // Want the exact intersection


      // Final Check (if there is overlap with our timestep)
      bool x_y_inter, y_z_inter, z_x_inter;
      x_y_inter = tx_min <= ty_max  && ty_min <= tx_max;
      y_z_inter = ty_min <= tz_max  && tz_min <= ty_max;
      z_x_inter = tz_min <= tx_max  && tx_min <= tz_max;

      // bool x_sat = tx_min <= time_step && time_step <= tx_max;
      // bool y_sat = ty_min <= time_step && time_step <= ty_max;
      // bool z_sat = tz_min <= time_step && time_step <= tz_max;

      // If we do not share a common interval, no collision ever
      if(!x_y_inter || !y_z_inter || !z_x_inter){
        return;
      }

      // find exactly when we hit our box
      double t_hitbox = std::max( std::max( tx_min, ty_min), tz_min );

      // If we hit it in this timestep
      if ( 0 <= t_hitbox && t_hitbox <= time_step){

        if(isLeaf()){

          tv.push_back(getTriangle());

        }else{

          // Try to iterate through my children
          right_volume->getTriangles(r, time_step, tv);
          left_volume->getTriangles(r, time_step, tv);

        }
      
      }//if
        
    }//end



    
  private:

    BVHNode * left_volume;    // All NULL if a leaf node
    BVHNode * right_volume;   
    BoundingBox bbox;         
    Triangle * tri_leaf; // NULL unless a leaf node
    
};

#endif
