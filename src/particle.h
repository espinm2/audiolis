#ifndef PARTICLE_H
#define PARTICLE_H


#include <cassert>
#include <iostream>
#include <cmath>
#include <string>

#include "vectors.h"

class Particle
{
    public:

      // defualt constructor
      Particle(){
        position     = Vec3f(0,0,0);
        oldPosition  = Vec3f(0,0,0);
        ampage       = 0;
        splits       = 0;
      }

      // constructor
      Particle(const Vec3f & pos, const Vec3f & old, const Vec3f c, double a, int s){
          position    = pos;
          oldPosition = old;
          ampage      = a;
          splits      = s;
          center      = c;
      }

      // Accessors
      const   Vec3f & getPos()        const { return position; }
      const   Vec3f & getOldPos()     const { return oldPosition;}
      const   Vec3f & getCenter()     const { return center; }
              Vec3f getDir()          const { Vec3f res = (position-center); res.Normalize();
                                              return res; }
      double  getAmp()               const { return ampage; }
      int     getSplit()             const { return splits; }


      // Modifiers
      void setPos     (const Vec3f & pos) { position = pos; }
      void setOldPos  (const Vec3f & pos) { oldPosition = pos; }
      void setCenter  (const Vec3f & pos) { center = pos; }
      void setAmp     (const double & a)  { ampage = a; }
      void setSplit   (const int    & s)  { splits = s; }

      // Debugging Functions
      friend std::ostream& operator<<(std::ostream &, const Particle &);

private:

      // Rep
      Vec3f position;
      Vec3f oldPosition;
      Vec3f center;
      double   ampage;
      int      splits;
};


// Printing particle inline
inline std::ostream & operator<<(std::ostream & leftOp, const Particle & rightOp){

    leftOp << "Loc " << &rightOp << " ";
    leftOp << "pos(" << rightOp.position.x()  << ", " << rightOp.position.y()  << ", " <<  rightOp.position.z()  << ") ";
    leftOp << "old(" << rightOp.oldPosition.x()  << ", " << rightOp.oldPosition.y()  << ", " <<  rightOp.oldPosition.z()  << ") ";
    leftOp << "dir(" << rightOp.getDir().x() << ", " << rightOp.getDir().y() << ", " <<  rightOp.getDir().z() << ") ";
    leftOp << "cen(" << rightOp.center.x()    << ", " << rightOp.center.y()    << ", " <<  rightOp.center.z()    << ") ";
    leftOp << "amp " << rightOp.ampage      << " ";
    leftOp << "spl " << rightOp.splits;
    return leftOp;
}

#endif // PARTICLE_H
