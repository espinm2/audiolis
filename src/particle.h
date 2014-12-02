#ifndef PARTICLE_H
#define PARTICLE_H

#include <glm/glm.hpp>
#include <cassert>
#include <iostream>
#include <string>
#include <cmath>
#include <string>

class Particle {

public:

  // Constructor //////////////////////////////////////////////////////////////
  Particle(){
    position     = (glm::vec3) NULL; 
    oldPosition  = (glm::vec3) NULL; 
    center       = (glm::vec3) NULL; 
    hitNorm      = (glm::vec3) NULL;

    wattage      = 0;
    freqency     = 0;
    splits       = 0;

    timeLeft     = 0;
    iterations   = 0;

    materialHit =  "";
  }

  Particle(const glm::vec3 & pos, const glm::vec3 & old, 
      const glm::vec3 cen, double amp, double freq, int s){

      // Remember you should set stepsLeft + iterations
      position    = pos;
      oldPosition = old;
      center      = cen;

      wattage      = amp;
      freqency    = freq;
      splits      = s;

      timeLeft    = 0;
      iterations  = 0;

      materialHit =  "";

  }

  // Accessors ////////////////////////////////////////////////////////////////
  const   glm::vec3 & getPos()        const { return position; }
  const   glm::vec3 & getOldPos()     const { return oldPosition;}
  const   glm::vec3 & getCenter()     const { return center; }
  const   glm::vec3 & getHitNorm()    const { return hitNorm; }

  double  getWatt()     const { return wattage; }
  double  getFreq()     const { return freqency; }
  int     getSplit()    const { return splits; }
  float   getTimeLeft() const { return timeLeft; }
  int     getIter()     const { return iterations; }

  glm::vec3 getDir() const { 

    // Computes the direction
    glm::vec3 res = (position-center); 
    res = glm::normalize(res);
    return res;

  }

  std::string     getMaterialHit(){ return materialHit; }

  // Modifiers ////////////////////////////////////////////////////////////////
  void setPos     (const glm::vec3 & pos) { position = pos; }
  void setOldPos  (const glm::vec3 & pos) { oldPosition = pos; }
  void setCenter  (const glm::vec3 & pos) { center = pos; }
  void setHitNorm (const glm::vec3 & pos) { hitNorm = pos; }

  void setWatt    (const double & a)  { wattage = a; }
  void setFreq    (const double & f)   {freqency = f;}

  void setSplit   (const int    & s)  { splits = s; }

  void setTime    (const float  & t)  { timeLeft = t ; }
  void decTime    (const float  & t)  { timeLeft = timeLeft - t; }
  void incIter    () { iterations++; }

  void setMaterial (const std::string & name ) { materialHit = name; } 

  // Debugging Functions //////////////////////////////////////////////////////
  friend std::ostream& operator<<(std::ostream &, const Particle &);

private:

  // Rep Geometery ////////////////////////////////////////////////////////////
  glm::vec3 position;       // Position to use for rendering
  glm::vec3 oldPosition;    // Position before rendering
  glm::vec3 center;         // Where the epi-center is
  glm::vec3 hitNorm;        // Angle in which we hit the mesh
  
  // Rep Sound ////////////////////////////////////////////////////////////////
  double    wattage;         // Used to calc Power of wave
  double    freqency;       // What freqency this particle represents
  int       splits;         // How many times our particle split

  float     timeLeft;       // How much time left until you hit wall
  int       iterations;     // How many iterations have I been around for

  std::string materialHit;  // Name of the material I will hit next

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

    leftOp << "amp " << rightOp.wattage << " ";

    leftOp << "spl " << rightOp.splits << " ";

    leftOp << "stp " << rightOp.timeLeft << " ";

    leftOp << "itr " << rightOp.iterations;

    return leftOp;
}

#endif // PARTICLE_H
