#include "mask.h"
#define RADIUS_PARTICLE_WAVE 0.1f

void Mask::renderOutline( std::vector<VBOPosNormalColor> & outline_verts){
  std::cout << "ERROR: NOT IMPLEMENTED OUTLINE" << std::endl;
  assert(false);

}

void Mask::renderCost( std::vector<VBOPosNormalColor> & cost_verts){

  // FIXME as of right now I am visualizing the distance to split
  
  assert(costVector.size() == maskParticles.size()); // Safety

  
  for( unsigned int i = 0; i < maskParticles.size(); i++ ){

    if(maskParticles[i] == NULL)
      continue;


    Particle * cur = maskParticles[i];
    int cost = costVector[i];
  
    double val= cost / (1.2 * RADIUS_PARTICLE_WAVE * 1000000); 


    glm::vec4 happyColor =  GiveHeapMapping(val);

    // Pushing a line segement from this point to center for happyness ////
    cost_verts.push_back(
        VBOPosNormalColor(
          maskCenter->getPos(), 
          maskCenter->getDir(), 
          happyColor));

     happyColor.a = 0;

    cost_verts.push_back(
        VBOPosNormalColor(
          cur->getPos(), 
          cur->getDir(), 
          happyColor));
  }//for
}//func



bool Mask::resSpit(std::vector<glm::vec3> & newPartPos){

  bool split_happened = false;

  for( int i = 0; i < maskParticles.size(); i++ ){
  
    Particle * curOuter = maskParticles[i];

    if(curOuter == NULL){
      continue;

    }
  
    float dist_center_outer = glm::distance(curOuter->getOldPos(), maskCenter->getOldPos());

    // if that edge is too streched out
    if( dist_center_outer > 2 * RADIUS_PARTICLE_WAVE  ){

      // Average both vectors postions
      glm::vec3 posA = maskCenter->getOldPos();
      glm::vec3 posB = curOuter->getOldPos();
      glm::vec3 posNew = ((float) 0.5 ) * (posA + posB);

      newPartPos.push_back(posNew);
      split_happened = true;
    }
      
  } // for each edge

  return split_happened;

  // Blow if working resSplit for projection
  /*
  bool split_happened = false;

  for( int i = 0; i < maskParticles.size(); i++ ){
  
    Particle * curOuter = maskParticles[i];

    if(curOuter == NULL){
      continue;

    }
  
    // if that edge is too streched out
    if(costVector[i] > RADIUS_PARTICLE_WAVE * 1000 ){


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
  */

}




void Mask::debugPrint(){

  std::cout << "============================\n";
  std::cout << "Center " << &(*maskCenter) << "\n";

  std::cout << "Particle Mask [";
  for(Particle * p : maskParticles){
    if( p != NULL)
      std::cout << &(*p) << "\t";
    else
      std::cout << "NULL    " << "\t";
  }
  std::cout << "]\n";

  std::cout << "Particle Cost [";
  for(int c : costVector){
    if( c == 1000000)
      std::cout << "INF      " << "\t";
    else
    std::cout << c << "\t\t";
  }

  std::cout << "]\n";

}
