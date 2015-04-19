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
  
    double val= cost / (RADIUS_PARTICLE_WAVE * 1000); 


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


double maskArea(){
  
  // Purpose:
  //   Returns the area that this mask encompasses this is used to find 
  //   out how many particles we expect within the space
  //  Args:
  //    None
  
  assert(false);
  return 0.0;

}

bool Mask::resSpit(std::vector<glm::vec3> & newPartPos){
  // Purpose:
  //  Returns true or false if the current space needs to split
  //  if true it will also add new particle positions to  newPartPos 
  //  Args:
  //    newPartPos : a vector where we keep points, passed by ref

  // TODO IMPLEMENT THIS FUNCTION TO ALLOW FOR SPITS 
  // THAT CAN SLIDE BETWEEN EACHOTHER !!!
  assert(false);
  
  
  // Code might be useful for projecting back on the scene
  // // Project this new particle back onto the sphere
  // float radius = glm::distance(maskCenter->getCenter(), maskCenter->getOldPos());
  // glm::vec3 newDir = glm::normalize(posNew - maskCenter->getCenter());
  // posNew = maskCenter->getCenter() + radius * newDir;

  return false;

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
