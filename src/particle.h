#ifndef PARTICLE_H
#define PARTICLE_H

#include <glm/glm.hpp>
#include <cassert>
#include <iostream>
#include <string>
#include <cmath>
#include <string>

typedef unsigned int uint;

class Particle {

public:

  // Constructor //////////////////////////////////////////////////////////////
  Particle(){
    position     = (glm::vec3) NULL; 
    oldPosition  = (glm::vec3) NULL; 
    center       = (glm::vec3) NULL; 
    direction    = (glm::vec3) NULL;

    wattage      = 0;
    freqency     = 0;
    splits       = 0;

    iterations   = 0;
    _isDead      = false;
    collisionSteps = 0; 

  }

  Particle(const glm::vec3 & pos, const glm::vec3 & old, 
      const glm::vec3 cen, double watts, double freq, int s){

    // Remember you should set stepsLeft + iterations
    position    = pos;
    oldPosition = old;
    center      = cen;
    wattage      = watts;
    freqency    = freq;
    splits      = s;
    iterations  = 0;
    _isDead      = false;
    collisionSteps = 0; 
    direction = glm::normalize( old - cen );
  }

  // Accessors /////////////////////n///////////////////////////////////////////
  const   glm::vec3 & getPos()        const { return position; }
  const   glm::vec3 & getOldPos()     const { return oldPosition;}
  const   glm::vec3 & getCenter()     const { return center; }
  const   glm::vec3 & getDir()        const { return direction; }

  double  getWatt()                   const { return wattage; }
  double  getFreq()                   const { return freqency; }

  int     getIter()                   const { return iterations; }
  int     getSplit()                  const { return splits; }

  bool    isDead()                    const { return _isDead; }
  bool    isAlive()                   const { return !_isDead; }

  uint    getCollisionSteps()            const { return collisionSteps; }

  // Modifiers ////////////////////////////////////////////////////////////////
  void setPos     (const glm::vec3 & pos) { position = pos; }
  void setOldPos  (const glm::vec3 & pos) { oldPosition = pos; }
  void setCenter  (const glm::vec3 & pos) { center = pos; }
  void setDir     (const glm::vec3 & pos) { direction = pos; }

  void setWatt    (const double & a)  { wattage = a; }
  void setFreq    (const double & f)   {freqency = f;}

  void setIter    (const int & i)      { iterations = i; }
  void incIter    () { iterations++; }
  void setSplit   (const int    & s)  { splits = s; }
  void kill       () { _isDead = true; }
  void setCollisionSteps( uint i){ assert( i >= 0 ); collisionSteps = i; }

  // Debugging Functions //////////////////////////////////////////////////////
  friend std::ostream& operator<<(std::ostream &, const Particle &);

private:

  // Rep Geometery ////////////////////////////////////////////////////////////
  glm::vec3 position;       // Position to use for rendering
  glm::vec3 oldPosition;    // Position before rendering
  glm::vec3 center;         // Where the epi-center is
  glm::vec3 direction;      // The directon we are traveling

  // In Simulation ////////////////////////////////////////////////////////////
  int       splits;         // How many times our particle split
  int       iterations;     // How many iterations have I been around for
  bool      _isDead;        // Replacement for the delete mask
  uint      collisionSteps; // Steps before we encounter collision

  // Rep Sound ////////////////////////////////////////////////////////////////
  double    wattage;         // Used to calc Power of wave
  double    freqency;       // What freqency this particle represents

};



// Printing particle inline
inline std::ostream & operator<<(std::ostream & leftOp, 
    const Particle & rightOp){

    leftOp << "Loc " << &rightOp << " ";

    leftOp << "isAlive " << rightOp.isAlive() << " ";

    leftOp << "pos(" << rightOp.position.x  << 
      ", " << rightOp.position.y  << ", " <<  rightOp.position.z  << ") ";

    leftOp << "old(" << rightOp.oldPosition.x  << 
      ", " << rightOp.oldPosition.y  << ", " <<  rightOp.oldPosition.z  << ") ";

    leftOp << "dir(" << rightOp.getDir().x << 
      ", " << rightOp.getDir().y << ", " <<  rightOp.getDir().z << ") ";

    leftOp << "cen(" << rightOp.center.x    << 
      ", " << rightOp.center.y    << ", " <<  rightOp.center.z    << ") ";

    leftOp << "amp " << rightOp.wattage << " ";

    leftOp << "spl " << rightOp.splits << " ";

    leftOp << "itr " << rightOp.iterations;

    return leftOp;
}

#endif // PARTICLE_H
