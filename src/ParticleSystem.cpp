#include "ParticleSystem.h"

#include <glm/glm.hpp>
#include <iostream>
#include <cmath>

#include "boundingbox.h"
#include "MersenneTwister.h"
#include "particle.h"
#include "argparser.h"
#include "ray.h"
#include "hit.h"
#include "triangle.h"
#include "edge.h"
#include "vertex.h"
#include "hash.h"
#include "mesh.h"
#include "collision_utils.h"


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
  for( int i = 0; i < 10000; i++){
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

    Particle * p = new Particle(pos,cursor,cursor,100,0,1000);
    setStepBeforeCollision(p);

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


  // Stuff for calc
  glm::vec3 oldPos = p->getOldPos();
  glm::vec3 dir = p->getDir();

  if( p->getSteps() > 0 ){
  
    glm::vec3 newPos( oldPos + dir * args->timestep );

    p->setPos(newPos);
    p->decSteps(); // count down
  
  }else{
  
  
  }

}

void ParticleSystem::moveCursor( const float & dx, 
    const float & dy, const float & dz ){

  // TODO FIX
  cursor+= glm::vec3(dx,dy,dz);

}

void ParticleSystem::setStepBeforeCollision(Particle * &p){
  // Assuming that directiona and center work
  
  // Create a ray
  Ray r(p->getPos(), p->getDir());
  Hit h;
  bool hitTriangle;
  bool backface = false;


  // For each triangle we will check which we hit
  for (triangleshashtype::iterator iter = mesh->triangles.begin();
       iter != mesh->triangles.end(); iter++) {

    Triangle *t = iter->second;
    glm::vec3 a = (*t)[0]->getPos();
    glm::vec3 b = (*t)[1]->getPos();
    glm::vec3 c = (*t)[2]->getPos();    

    if(triangle_intersect(r,h,a,b,c,backface))
      hitTriangle = true;
  }

  if(hitTriangle){
    float time =  h.getT();
    p->setSteps( (int)(time / args->timestep) );

  }else{
    p->setSteps( 1000000 );

  }
}

