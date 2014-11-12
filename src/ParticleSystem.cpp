#include "ParticleSystem.h"

#include <glm/glm.hpp>
#include <iostream>
#include <cmath>

#include "boundingbox.h"
#include "MersenneTwister.h"
#include "particle.h"
#include "argparser.h"

// Used for update
typedef std::vector<Particle *>::iterator ParticleIter;

ParticleSystem::~ParticleSystem(){
  // Just delete all the particles we made
  for(int i = 0; i < particles.size(); i++)
    delete particles[i];
}


void ParticleSystem::load(){
  
  // Initaite the cursor 
  glm::vec3 centerScene;

  bbox->getCenter(centerScene);

  cursor = glm::vec3(centerScene.x, centerScene.y, centerScene.z);

  // Using cursor
  MTRand randomGen;
  double s = 0.05;

  // Create Box of ranodm points
  for( int i = 0; i < 1000; i++){
    // Find x,y,z
    float x = cursor.x - s/2.0 + (float) randomGen.rand(s);
    float y = cursor.y - s/2.0 + (float) randomGen.rand(s);
    float z = cursor.z - s/2.0 + (float) randomGen.rand(s);

    glm::vec3 pos(x,y,z);

    // Project into a circle
    float radius = s * sqrt(2.0);
    
    glm::vec3 dir = pos - cursor;

    dir = glm::normalize(dir);

    pos = cursor + dir * radius;

    Particle * p = new Particle(pos,cursor,cursor,100,0);

    // put particle there
    particles.push_back(p);

  }

}

void ParticleSystem::update(){

  // TODO Implement collision detection + splits

  for(ParticleIter iter = particles.begin(); iter != particles.end(); iter++){

    Particle * curPart = (*iter);

    // New position is now old
    curPart->setOldPos(curPart->getPos());

    moveParticle(curPart);
  
  }

}

void ParticleSystem::moveParticle(Particle * &p){

  glm::vec3 oldPos = p->getOldPos();
  glm::vec3 dir = p->getDir();


  // TODO Figure out if we want to use center for time, or baby ray tech
  glm::vec3 newPos( oldPos + dir * args->timestep );

  p->setPos(newPos);

}

void ParticleSystem::moveCursor( const float & dx, 
    const float & dy, const float & dz ){

  // TODO FIX
  cursor+= glm::vec3(dx,dy,dz);

}
