#include "mask.h"

void Mask::renderOutline( std::vector<VBOPosNormalColor> & outline_verts){
  std::cout << "ERROR: NOT IMPLEMENTED OUTLINE" << std::endl;
  assert(false);

}

void Mask::renderCost( std::vector<VBOPosNormalColor> & cost_verts){

  assert(costVector.size() == maskParticles.size()); // Safety

  if( size_of_mask < 6 )
    return;
  
  for( unsigned int i = 0; i < maskParticles.size(); i++ ){

    Particle * cur = maskParticles[i];
    int cost = costVector[i];
  
    double val= cost / (1.6*1000* 0.1); // <--------------------------------------- REALLY BAD CODE :( FORGIVE ME
    glm::vec4 happyColor =  GiveHeapMapping(val);

    // Pushing a line segement from this point to center for happyness ////
    cost_verts.push_back(VBOPosNormalColor(maskCenter->getPos(), maskCenter->getDir(), happyColor));
    cost_verts.push_back(VBOPosNormalColor(cur->getPos(), cur->getDir(), happyColor));
  }//for
}//func
