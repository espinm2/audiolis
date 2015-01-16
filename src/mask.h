#ifndef MASK_H
#define MASK_H

#include "particle.h"
#include "render_utils.h"
#include  <vector>


// This is a tempurary datastructure that will create
class Mask{

  public:

    // Constructor
    Mask(Particle * cen, std::vector <Particle *> maskP){
      maskCenter = cen;
      maskParticles = maskP;
    }
    
    // Accessors
    const std::vector<Particle *> & getMaskParticles() { return maskParticles;}
    const Particle * getCenter() { return maskCenter;}

    void renderOutline( std::vector<VBOPosNormalColor> & outline_verts);
    void renderCost( std::vector<VBOPosNormalColor> & cost_verts);


  private:
    Particle * maskCenter;
    std::vector< Particle *> maskParticles;

};



#endif // MASK_H
