#include "ParticleSystem.h"

#include <glm/glm.hpp>
#include <iostream>
#include <cmath>
#include <limits>

#include "boundingbox.h"
#include "geometry_utils.h"
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
#include "render_utils.h"


// FIXME ////////////////////////////////////////////
// [ ] Particle waves should happen with distance
// [ ] Get low resolution of a mesh to use (tri count)
// [ ] when you split, fix behind timestep issue
// [ ] Change the radius_particle_wave depending how many spits done
// [ ] ReIntroduce splits

#define MAX_ITERATIONS        600   // When to kill particles
#define RADIUS_PARTICLE_WAVE  0.01 // Radius of hex shape
#define RADIUS_INIT_SPHERE    0.01  // Radius of source sphere
#define NUM_INIT_PARTICLES    1 // Inital Number of Particle
#define INITAL_AMPLATUDE      100   // Amp we start off with
#define SPLIT_AMOUNT          6     // What sized polygon we split into
#define MIN_AMP               10    // When particles should die

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
  /*
   * Input : None
   * Output: This function will update our particle simuations
   * Asumpt: There are particles to move
   * SideEf: Updates postition of particles/ removes particles
   */

  // Hold new particles from split
  std::vector<Particle *> newParticles;

  // Marked for removal mask, 1 == delete, 0 == keep
  std::vector<int>deleteMask (particles.size(), 0);
  unsigned int maskIndex = 0;

  for(ParticleIter iter = particles.begin(); iter != particles.end();){

     // This current Particle
     Particle * curPart = (*iter);
     curPart->setOldPos(curPart->getPos());


    // Are we below a threhold just kill and move to another
    if(curPart->getAmp() < MIN_AMP){

      // Kill this partcile and move to next
      deleteMask[maskIndex++] = 1;
      iter++;
      continue;

    }

    // Particles are beyond a threshold init a split
    if(shouldSplit(curPart)){

        // Get new particles to be made
        std::vector<Particle *> splitParticles;
        particleSplit(curPart, splitParticles);


        // Move them a timestep + add to new list
        for(int i = 0; i < splitParticles.size(); i++){
          moveParticle(splitParticles[i], args->timestep);
          newParticles.push_back(splitParticles[i]);
        }

      deleteMask[maskIndex++] = 0; // kills center
      iter++;


    }else{
    
      // Update postiton and move to next particle
      moveParticle(curPart,args->timestep);


      deleteMask[maskIndex++] = 0;
      iter++;
    
    }


  } //forloop


  // Deletetion step
  for( unsigned int i = 0 ; i < particles.size(); i++){
      // Keep if 0, else delete
      if(deleteMask[i] == 1){

          if(!newParticles.empty()){

              // Put in new particle to fill gap
              delete particles[i];
              particles[i] = newParticles.back();
              newParticles.pop_back();

          }else{

              // there is stuff to push off
              if(i != particles.size()-1){

                // Pop off back of vector to fill the gap
                delete particles[i];
                particles[i] = particles.back();
                particles.pop_back();

              }else{

                // Just delete the last element, nothing need be poped
                delete particles[i];
                particles.pop_back();

              }
          }
      }
  }

  // Add into the main vector those new particles yet added
  for( unsigned int i = 0; i < newParticles.size(); i++)
      particles.push_back(newParticles[i]);
}

bool ParticleSystem::moveParticle(Particle * p, float timestep){
  /*
   * Input : Particle ptr
   * Output: That particle moved a timestep
   * Asumpt: That particle exisit
   * SideEf: Changes p->position + might chug  when we change dir
   * SideEf: Changes center of particles
   *         
   */

  // Stuff for calc
  float time_until_impact = p->getTimeLeft();
  float time_after_impact = timestep - p->getTimeLeft();

  glm::vec3 oldPos = p->getOldPos();
  glm::vec3 dir = p->getDir();

  // We didn't hit an object in this interval of time
  if(time_until_impact - timestep > 0){

    glm::vec3 newPos( oldPos + dir * timestep );
    p->setPos(newPos);
    p->decTime(timestep);
  
  }else{

    // Get  Times
    float time_until_impact = p->getTimeLeft();
    float time_after_impact = timestep - p->getTimeLeft();

    // We where we hit in space
    glm::vec3 impactPos(oldPos + dir * time_until_impact);

    // Get the new center to change direction
    glm::vec3 mir_dir = MirrorDirection(p->getHitNorm(), p->getDir());
    mir_dir = mir_dir * (float)(-1.0);
    float radius = glm::distance(p->getCenter(), impactPos);
    
    // Move up to line
    p->setCenter(impactPos + mir_dir * radius);
    p->setOldPos(impactPos);
    calcMeshCollision(p); // new timeLeft Case a) we move a little up

    moveParticle(p, time_after_impact); // this dude will move
  }
  // up inter
  p->incIter();
}

void ParticleSystem::moveCursor( const float & dx, 
    const float & dy, const float & dz ){

  cursor+= glm::vec3(dx,dy,dz);

}

void ParticleSystem::splitAllParticles(){

  // Dead to me!
  // for(int i = 0; i < particles.size(); i++){
  //   particleSplit(particles[i]);
  // }

  // for(int i = 0; i < newParticles.size(); i++){
  //   particles.push_back(newParticles[i]);
  // }

  // newParticles.clear();

}

void ParticleSystem::particleSplit(Particle * &p,
 std::vector<Particle *> &vec){

  // Side effects: fills vec with new particles
  // Note: Particles are not updated at this step

    // Where I will store new particles
    std::vector< glm::vec3> newPart;

    // Get hex shape on plane
    circle_points_on_plane(p->getOldPos(), 
        p->getDir(), RADIUS_PARTICLE_WAVE, SPLIT_AMOUNT, newPart);
 
    // Project back on sphere
    cirlce_point_on_sphere(p->getCenter(),
        glm::distance( p->getOldPos(), p->getCenter()),newPart);

    // For each calculated pos, make particle
    for(int i = 0; i < newPart.size(); i++){
    
      Particle * s = new Particle(newPart[i],newPart[i],
        p->getCenter(),INITAL_AMPLATUDE, p->getSplit() + 1);
      calcMeshCollision(s);
    
      // put particle there to be "moved" when its their turn
      vec.push_back(s);
    }// for

}

void ParticleSystem::calcMeshCollision(Particle * &p){
  // Input  :  Give a particle
  // Output :  None
  // Assumpt:  That direction and center are setup
  // SideEff:  Sets timeLeft until collision
  
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

  // Get the closest  hit
  if(hitTriangle){
    p->setTime(h.getT());
    p->setHitNorm(h.getNormal());
  }else{
    // We don't hit any triangles
    p->setTime(10000000); 
  }
}


void ParticleSystem::createInitWave(){

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
  
    Particle * p = new Particle(pos,cursor,cursor,INITAL_AMPLATUDE,0);
    calcMeshCollision(p);
  
    // put particle there
    particles.push_back(p);
  
  }

}

bool ParticleSystem::shouldSplit(Particle * &p){

  // tag
  glm::vec3 pos = p->getOldPos();
  float nearestDistance = 100000;
  float threshold = RADIUS_PARTICLE_WAVE * 2;

  // for each particle in the system
  for(int i = 0; i < particles.size(); i++){

    if(particles[i] == p)
      continue;

    float dist = glm::distance(p->getOldPos(), particles[i]->getOldPos());
    if(nearestDistance > dist)
      nearestDistance = dist;
  }

  // I now have nearest distance
  return( fabs( nearestDistance - threshold) <= EPSILON || 
      nearestDistance > threshold);
}

