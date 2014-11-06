#include "ParticleSystem.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "vectors.h"
#include "boundingbox.h"

ParticleSystem::~ParticleSystem(){

  //TODO Implement

}


void ParticleSystem::load(){
  
  // Initaite the cursor 
  glm::vec3 centerScene;
  bbox->getCenter(centerScene);
  cursor = Vec3f(centerScene.x, centerScene.y, centerScene.z);

}

void ParticleSystem::update(){
  
  //TODO Implement

}

void ParticleSystem::moveCursor( const double & dx, const double & dy, const double & dz ){

  cursor = Vec3f(cursor.x() + dx,
                 cursor.y() + dy,
                 cursor.z() + dz);

}
