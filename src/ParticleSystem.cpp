#include "ParticleSystem.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <cmath>
#include "vectors.h"
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
  cursor = Vec3f(centerScene.x, centerScene.y, centerScene.z);

  // Using cursor
  MTRand randomGen;
  double s = 0.05;

  // Create Box of ranodm points
  for( int i = 0; i < 1000; i++){
    // Find x,y,z
    double x = cursor.x() - s/2.0 + randomGen.rand(s);
    double y = cursor.y() - s/2.0 + randomGen.rand(s);
    double z = cursor.z() - s/2.0 + randomGen.rand(s);

    Vec3f pos(x,y,z);

    // Project into a circle
    double radius = s * sqrt(2.0);
    
    Vec3f dir = pos - cursor;
    dir.Normalize();

    pos = cursor + dir * radius;

    Particle * p = new Particle(pos,cursor,cursor,100,0);

    std::cout << dir.x() <<  " " << dir.y() << " " << dir.z() << std::endl;

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
    // curPart->setPos(Vec3f(0,0,0)); // New is cleared
    moveParticle(curPart);
  
  }

}

void ParticleSystem::moveParticle(Particle * &p){

  Vec3f oldPos = p->getOldPos();
  Vec3f dir = p->getDir();



  // TODO Figure out if we want to use center for time, or baby ray tech
  Vec3f newPos( oldPos + dir * args->timestep );

  p->setPos(newPos);

}

void ParticleSystem::moveCursor( const double & dx, 
    const double & dy, const double & dz ){

  cursor = Vec3f(cursor.x() + dx,
                 cursor.y() + dy,
                 cursor.z() + dz);
}
