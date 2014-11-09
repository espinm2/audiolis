#include "ParticleSystem.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include "vectors.h"
#include "boundingbox.h"
#include "MersenneTwister.h"
#include "particle.h"

ParticleSystem::~ParticleSystem(){
  // Just delete all the particles we made
  for(int i = 0; i < particles.size(); i++)
    delete particles[i];
}


void ParticleSystem::load(){
  
  // Initaite the cursor 
  glm::vec3 centerScene;
  bbox->getCenter(centerScene);
  cursor = Vec3f(centerScene.x, centerScene.y, centerScene.z);

  // Using cursor
  MTRand randomGen;
  double s = 1;

  // Create Box of ranodm points
  for( int i = 0; i < 1000; i++){
    // Find x,y,z
    double x = cursor.x() - s/2.0 + randomGen.rand(s);
    double y = cursor.y() - s/2.0 + randomGen.rand(s);
    double z = cursor.z() - s/2.0 + randomGen.rand(s);

    Vec3f pos(x,y,z);
    Particle * p = new Particle(pos,cursor,cursor,100,0);
    std::cout << *p << std::endl;

    // put particle there
    particles.push_back(p);

  }

}

void ParticleSystem::update(){
  
  //TODO Implement

}

void ParticleSystem::moveCursor( const double & dx, 
    const double & dy, const double & dz ){

  cursor = Vec3f(cursor.x() + dx,
                 cursor.y() + dy,
                 cursor.z() + dz);
}
