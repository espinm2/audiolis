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
#include "render_utils.h"
#include "split_utils.h"


// FIXME ////////////////////////////////////////////
// [ ] Particle waves should happen with distance
// [ ] Get low resolution of a mesh to use (tri count)
// [ ] when you split, fix behind timestep issue
// [ ] Change the radius_particle_wave depending how many spits done
// [ ] ReIntroduce splits

#define MAX_ITERATIONS 600 // refactor with below
#define RADIUS_PARTICLE_WAVE 0.001 //TODO hex created size should very
#define RADIUS_INIT_SPHERE 0.01
#define NUM_INIT_PARTICLES  100
#define NUM_INIT_LIFESPACE  1000 // refactor
#define SPLIT_AMOUNT 6

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
}

void ParticleSystem::update(){


  // for(ParticleIter iter = particles.begin(); iter != particles.end();){
  for(int i = 0; i < particles.size();){

    //Particle * curPart = (*iter);
    Particle * curPart = particles[i];


    if(moveParticle(curPart)){
      i++;
    }else{

        //TODO  Make this erease more effiecent
        Particle * backParticle = particles.back();

        if( curPart == backParticle ){

          delete curPart; // kill what p is pointing to
          particles.pop_back();


        }else{


          delete curPart;

          particles[i] = backParticle;

          particles.pop_back();


      }
    }
  }

}

bool ParticleSystem::moveParticle(Particle * p){


  if(p->getIter() > MAX_ITERATIONS){


    return false;
  }

  // up inter
  p->incIter();
  
  // New position is now old
  p->setOldPos(p->getPos());

  // Stuff for calc
  glm::vec3 oldPos = p->getOldPos();
  glm::vec3 dir = p->getDir();

  if( p->getSteps() > 0 ){
  

    glm::vec3 newPos( oldPos + dir * args->timestep );

    p->setPos(newPos);
    p->decSteps(); // count down

    
  
  }else{

    dir = MirrorDirection(p->getHitNorm(), p->getDir());
    dir = dir * (float)(-1.0);

    float radius = glm::distance(p->getCenter(), p->getOldPos());

    p->setCenter(oldPos + dir * radius);

    setStepBeforeCollision(p);
    // New upated for new center/dir
    glm::vec3 newdir = p->getDir();
    glm::vec3 newPos( oldPos + newdir * args->timestep );
    p->setPos(newPos);
    p->decSteps(); // count down
  
  }

  return true;

}

void ParticleSystem::moveCursor( const float & dx, 
    const float & dy, const float & dz ){

  cursor+= glm::vec3(dx,dy,dz);

}

void ParticleSystem::splitAllParticles(){

  
  for(int i = 0; i < particles.size(); i++){
    particleSplit(particles[i]);
  }

  for(int i = 0; i < newParticles.size(); i++){
    particles.push_back(newParticles[i]);
  }

  newParticles.clear();

}

//tag
void ParticleSystem::particleSplit(Particle * &p){
  // Side effects: Adds new particles to end of particles list

    // Where I will store new particles
    std::vector< glm::vec3> newPart;

    circle_points_on_plane(p->getOldPos(), 
        p->getDir(), RADIUS_PARTICLE_WAVE, SPLIT_AMOUNT, newPart);
 
    cirlce_point_on_sphere(p->getCenter(),
        glm::distance( p->getOldPos(), p->getCenter()),newPart);

    for(int i = 0; i < newPart.size(); i++){
    
      // FIXME this could cause a problem in future
      // If you ever use oldPos or newPos for any visual
      Particle * s = new Particle(newPart[i],newPart[i],
          p->getCenter(),100,0, NUM_INIT_LIFESPACE);

      setStepBeforeCollision(s);
    
      // put particle there to be "moved" when its their turn
      newParticles.push_back(s);
    
    }// for
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
    p->setHitNorm(h.getNormal());

  }else{
    p->setSteps( 1000000 );

  }
}


void ParticleSystem::createParticleWave(){


  // Testing function to create circle in 3d space

  // Phonon Mapping to create particle wave
  // Using cursor
  MTRand randomGen;
  double s = RADIUS_INIT_SPHERE;
  
  // Create Box of ranodm points
  for( int i = 0; i < NUM_INIT_PARTICLES; i++){
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
  

    Particle * p = new Particle(pos,cursor,cursor,100,0,NUM_INIT_LIFESPACE);
    setStepBeforeCollision(p);
  
    // put particle there
    particles.push_back(p);
  
  }

}



