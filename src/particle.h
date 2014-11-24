#ifndef PARTICLE_H
#define PARTICLE_H

#include <glm/glm.hpp>
#include <cassert>
#include <iostream>
#include <cmath>
#include <string>

class Particle {

public:

  // defualt constructor
  Particle(){
    position     = glm::vec3(0,0,0);
    oldPosition  = glm::vec3(0,0,0);
    center       = glm::vec3(0,0,0);
    hitNorm      = (glm::vec3) NULL;
    ampage       = 0;
    splits       = 0;
    timeLeft     = 0;
    iterations   = 0;
  }

  // constructor
  Particle(const glm::vec3 & pos, const glm::vec3 & old, 
      const glm::vec3 c, double a, int s){

      // Remember you should set stepsLeft + iterations
      position    = pos;
      oldPosition = old;
      center      = c;
      ampage      = a;
      splits      = s;
      timeLeft    = 0;
      iterations  = 0;

  }

  // Accessors
  const   glm::vec3 & getPos()        const { return position; }
  const   glm::vec3 & getOldPos()     const { return oldPosition;}
  const   glm::vec3 & getCenter()     const { return center; }
  const   glm::vec3 & getHitNorm()    const { return hitNorm; }

  double  getAmp()               const { return ampage; }
  int     getSplit()             const { return splits; }
  float   getTimeLeft()          const { return timeLeft; }
  int     getIter()              const { return iterations; }

  glm::vec3 getDir() const { 
    glm::vec3 res = (position-center); 
    res = glm::normalize(res);
    return res;
  }


  // Modifiers
  void setPos     (const glm::vec3 & pos) { position = pos; }
  void setOldPos  (const glm::vec3 & pos) { oldPosition = pos; }
  void setCenter  (const glm::vec3 & pos) { center = pos; }
  void setHitNorm (const glm::vec3 & pos) { hitNorm = pos; }
  void setAmp     (const double & a)  { ampage = a; }
  void setSplit   (const int    & s)  { splits = s; }
  void setTime    (const float  & c)  { timeLeft = c;   }
  void decTime    (const double & t)  { timeLeft -= (float)t; }
  void incIter    () { iterations++; }

  // Debugging Functions
  friend std::ostream& operator<<(std::ostream &, const Particle &);

private:

  // Rep
  glm::vec3 position;       // Position to use for rendering
  glm::vec3 oldPosition;    // Position before rendering
  glm::vec3 center;         // Where the epi-center is
  glm::vec3 hitNorm;        // Angle in which we hit the mesh
  
  double    ampage;         // Used to calc Power of wave
  int       splits;         // How many times our particle split
  float     timeLeft;       // How much time left until you hit wall
  int       iterations;     // How many iterations have I been around for
};


// Printing particle inline
inline std::ostream & operator<<(std::ostream & leftOp, 
    const Particle & rightOp){

    leftOp << "Loc " << &rightOp << " ";

    leftOp << "pos(" << rightOp.position.x  << 
      ", " << rightOp.position.y  << ", " <<  rightOp.position.z  << ") ";

    leftOp << "old(" << rightOp.oldPosition.x  << 
      ", " << rightOp.oldPosition.y  << ", " <<  rightOp.oldPosition.z  << ") ";

    leftOp << "dir(" << rightOp.getDir().x << 
      ", " << rightOp.getDir().y << ", " <<  rightOp.getDir().z << ") ";

    leftOp << "cen(" << rightOp.center.x    << 
      ", " << rightOp.center.y    << ", " <<  rightOp.center.z    << ") ";

    leftOp << "hitNorm(" << rightOp.hitNorm.x << 
      ", " << rightOp.hitNorm.y << ", " << rightOp.hitNorm.z << ") ";

    leftOp << "amp " << rightOp.ampage << " ";

    leftOp << "spl " << rightOp.splits << " ";

    leftOp << "stp " << rightOp.timeLeft << " ";

    leftOp << "itr " << rightOp.iterations;

    return leftOp;
}

#endif // PARTICLE_H
