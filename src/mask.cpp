#include "mask.h"
#define RADIUS_PARTICLE_WAVE 0.1f

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
  
    double val= cost / (1.6*1000* RADIUS_PARTICLE_WAVE); // <--------------------------------------- REALLY BAD CODE :( FORGIVE ME
    glm::vec4 happyColor =  GiveHeapMapping(val);

    // Pushing a line segement from this point to center for happyness ////
    cost_verts.push_back(VBOPosNormalColor(maskCenter->getPos(), maskCenter->getDir(), happyColor));
    cost_verts.push_back(VBOPosNormalColor(cur->getPos(), cur->getDir(), happyColor));
  }//for
}//func



bool Mask::resSpit(std::vector<glm::vec3> & newPartPos){


  bool split_happened = false;

  for( int i = 0; i < maskParticles.size(); i++ ){
  
    Particle * curOuter = maskParticles[i];

    if(curOuter == NULL){
      continue;

    }
  
    // if that edge is too streched out
    if(costVector[i] > RADIUS_PARTICLE_WAVE * 1000 ){

      std::cout << "Edge to far, split init" << std::endl;

      // Average both vectors postions
      glm::vec3 posA = maskCenter->getOldPos();
      glm::vec3 posB = curOuter->getOldPos();
      glm::vec3 posNew = ((float) 0.5 ) * (posA + posB);
      glm::vec3 center = maskCenter->getCenter();

      // Get direction from the center
      glm::vec3 dirNew = glm::normalize(posNew - center);
    
      // Project on the imaginary sphere
      float dist = glm::distance(posA, center);
      posNew = center + dirNew * dist;

      newPartPos.push_back(posNew);
      split_happened = true;
    }
      
  } // for each edge

  return split_happened;

}
