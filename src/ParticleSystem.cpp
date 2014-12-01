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
#define INITAL_AMPLATUDE      10000000   // Amp we start off with
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

    // Particles are beyond a threshold init a split
    if(shouldSplit(curPart)){

        if(curPart->getAmp() < MIN_AMP){

          // Kill this partcile and move to next
          deleteMask[maskIndex++] = 1;
          iter++;
          continue;

        }

        // Get new particles to be made
        std::vector<Particle *> splitParticles;
        particleSplit(curPart, splitParticles);

        // change cur particle amp
        curPart->setAmp(curPart->getAmp() / (double)(SPLIT_AMOUNT + 1.0));

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

void ParticleSystem::particleSplit(Particle * &p,
 std::vector<Particle *> &vec){
  // Side effects: fills vec with new particles
  // Note: Particles are not updated at this step

    // Where I will store new particles
    std::vector< glm::vec3> newPart;
    
    // Get hex shape on plane
    circle_points_on_plane(p->getOldPos(), 
        p->getDir(), RADIUS_PARTICLE_WAVE, SPLIT_AMOUNT, newPart,args);

 
    // Project back on sphere
    cirlce_point_on_sphere(p->getCenter(),
        glm::distance( p->getOldPos(), p->getCenter()),newPart);

    // For each calculated pos, make particle
    for(int i = 0; i < newPart.size(); i++){
    
      Particle * s = new Particle(
          newPart[i], newPart[i], p->getCenter(),
          p->getAmp() / (double)(SPLIT_AMOUNT + 1.0), 
          p->getSplit() + 1);
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

  // p->setTime(10000000);  // <------------------------------REMOVE ME 
}


void ParticleSystem::createInitWave(){
  // Testing function to create circle in 3d space

  // Phonon Mapping to create particle wave
  // Using cursor
  double s = RADIUS_INIT_SPHERE;
  
  // Create Box of ranodm points
  for( int i = 0; i < NUM_INIT_PARTICLES; i++){
    // Find x,y,z
    float x = cursor.x - s/2.0 + (float) args->randomGen.rand(s);
    float y = cursor.y - s/2.0 + (float) args->randomGen.rand(s);
    float z = cursor.z - s/2.0 + (float) args->randomGen.rand(s);
  
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

  // DEBUG //
  return false;
  //  // tag
  //  glm::vec3 pos = p->getOldPos();
  //  float nearestDistance = 100000;
  //  float threshold = 3 * RADIUS_PARTICLE_WAVE;

  //  // for each particle in the system
  //  for(int i = 0; i < particles.size(); i++){

  //    if(particles[i] == p)
  //      continue;

  //    float dist = glm::distance(p->getOldPos(), particles[i]->getOldPos());
  //    if(nearestDistance > dist)
  //      nearestDistance = dist;
  //  }

  //  // I now have nearest distance
  //  return( fabs( nearestDistance - threshold) <= EPSILON || 
  //      nearestDistance > threshold);
}

void ParticleSystem::particleMerge(const Particle * &a, 
    const Particle * &b, Particle * &c){ 
  //TODO implement this function ( Solve this problem )
}


double ParticleSystem::absorbFunc(const std::string & materialName, 
    const double freq){

  // Insurance
  if( freq < 20 ){
    std::cout << "Freq too low" << std::endl;
    assert(false);
  }

  // This method will have the manually loaded absorb methods taken from
  // the phonon tracing paper imported into it!

  if( materialName == "concrete_wall"){

    if(4063 <= freq)
      return 0.21;

    else if( 2031 <= freq )
      return 0.00009686 * (freq - 2031) + 0.17;

    else if( 1015 <= freq )
      return 0.00000984 * ( freq - 1015 ) + 0.16;

    else if( 507 <= freq )
      return -0.00077165 * (freq - 507) +0.25;

    else if( 253 <= freq )
      return 0.00039708 * (freq - 253 )  + 0.15;

    else if( 126 <= freq )
      return 0.00070866 * (freq - 126 ) + 0.06;

    else 
      return 0.06;

  }else if(materialName == "plaster_ceiling"){

    if(1015 <= freq)
      return 0.05;

    else if(507 <= freq )
      return -0.00009843 * (freq - 507) + 0.1;

    else if(235 <= freq )
      return -0.00036765 * (freq - 235) + 0.2;

    else
      return 0.2;

  }else if(materialName == "ceramic_wall"){

    if(freq <= 2031)
      return 0.07;

    else if( 507 <= freq )
      return 0.00001312 * ( freq - 507 ) + 0.05;

    else if( 126 <= freq )
      return 0.00007874 * (freq - 126) + 0.02;

    else
      return 0.02;

  }else if(materialName == "bricked_wall"){

    if(4063 <= freq )
      return 0.03;

    else if(2031 <= freq )
      return 0.00000492 * (freq - 2031) + 0.02;

    else if(507 <= freq )
      return 0.02;
        
    else if(253 <= freq )
      return 0.00003937 * (freq - 253) + 0.01;

    else
      return 0.01;

  }else if(materialName == "double_window"){

    if(1015 <= freq)
      return 0.02;

    else if( 253 <= freq)
      return -0.00002625 * (freq - 253) + 0.04;

    else if (126 <= freq)
      return -0.00047244 * (freq - 126) + 0.1;

    else
      return 0.1;

  }else if(materialName == "pvc_floor"){

    if( 2031 <= freq )
      return 0.05;

    else if(507 <= freq )
      return 0.00002625(freq - 507) + 0.01;

    else
      return 0.1;

  }else if(materialName == "carpated_floor"){
    
    if( 4063 <= freq)
      return 0.35;

    else if( 2031 <= freq )
      return 0.00007874 * (freq - 2031 ) + 0.19;

    else if( 507 <= freq )
      return 0.00009186 * (freq - 507) + 0.05;

    else
      return 0.00004107 * (freq-20) + 0.03;

  }else if(materialName == "absorber_parete"){

    if( 253 <= freq )
      return 0.9;

    else if(126 <= freq)
      return 0.00157480 * (freq - 126) + 0.7;

    else if(63.4 <= freq)
      return 0.00878594 * (freq - 63.4) + 0.15;

    else if(31.7 <= freq )
      return 0.00315457 * ( freq - 31.7 ) + 0.05;

    else
      return 0.05;
  
  }else{
    std::cout<< "Error, no correct material associated with an object\n";
    assert(false);
  }

}

