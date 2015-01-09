#include "ParticleSystem.h"

#include <glm/glm.hpp>
#include <iostream>
#include <cmath>
#include <limits>
#include <float.h>
#include <cstdlib>

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
#include "munkers.h"
#include <algorithm>


#define EPSILON 0.0001

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

  // Initiate bounding box
  bbox->getCenter(centerScene);

  // put  cursor in center of scene
  cursor = glm::vec3(centerScene.x, centerScene.y, centerScene.z);

  // Initalize simulation variables
  TIME_STEP             = args->timestep; // imported from args file
  VELOCITY_OF_MEDIUM    = 340; //implement to sound speed in m/s

  // simuation var init
  RADIUS_INIT_SPHERE    = 0.1; // Doesnt get changed to often
  NUM_INIT_PARTICLES    = args->num_init_particles; // imported from args
  MIN_WATTAGE           = 0.000000000002; // Doesnt change often
  MAX_ITERATIONS        = 6000; // Doesnt change often

  // split var init
  RADIUS_PARTICLE_WAVE  = 0.01; // Will be later removed for better split
  SPLIT_AMOUNT          = 6; // Will be later removed for better split


  // Debug function used to test things upon creations
  // debug();

}

void ParticleSystem::debug(){
  // Currently just testing if our merge function works as expected

  /*
  // Create particle template
  Particle * a = new Particle(glm::vec3(0,0,0), // NewPos
                              glm::vec3(0,0,0), // OldPos
                              glm::vec3(0,0,0), // Cen
                              0, // Watts
                              0, // Freq
                              0); // spits
                              a->setIter(0); // iterations
  */


  // Testing oldPos
  Particle * a = new Particle(glm::vec3(0,0,0), // NewPos
                              glm::vec3(0,0,0), // OldPos
                              glm::vec3(0,0,0), // Cen
                              0, // Watts
                              0, // Freq
                              0); // spits
                              a->setIter(0); // iterations

  Particle * b = new Particle(glm::vec3(0,0,0), // NewPos
                              glm::vec3(1,1,1), // OldPos
                              glm::vec3(0,0,0), // Cen
                              0, // Watts
                              0, // Freq
                              0); // spits
                              b->setIter(0); // iterations

  Particle * c = particlePairMerge(a,b);

  std::cout << &c << std::endl;


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
    if(curPart->getWatt() < MIN_WATTAGE ){ 

      // Kill this partcile and move to next
      deleteMask[maskIndex++] = 1;
      iter++;
      continue;

    }

    // Particles are beyond a threshold init a split
    if(particleSplitCheckAndMerger(curPart, deleteMask)){

        // should we even bother? are they just going to flitter out
        if(curPart->getWatt()/(SPLIT_AMOUNT+1.0) < MIN_WATTAGE ){
          deleteMask[maskIndex++] = 1;
          iter++;
          continue;
        }

        // Get new particles to be made
        std::vector<Particle *> splitParticles;
        particleSplit(curPart, splitParticles);

        // change cur particle amp
        curPart->setWatt(curPart->getWatt() / (double)(SPLIT_AMOUNT + 1.0));

        // Move them a timestep + add to new list
        for(int i = 0; i < splitParticles.size(); i++){
          moveParticle(splitParticles[i], TIME_STEP);
          newParticles.push_back(splitParticles[i]);
        }

      deleteMask[maskIndex++] = 0; // 1: kills center, 0: leave it alive

      iter++;


    }else{
    
      // Update postiton and move to next particle
      moveParticle(curPart,TIME_STEP);
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

bool ParticleSystem::moveParticle(Particle * p, double timestep){
  /*
   * Input : Particle ptr
   * Output: That particle moved a timestep
   * Asumpt: That particle exisit
   * SideEf: Changes p->position + might chug  when we change dir
   * SideEf: Changes center of particles
   *         
   */

  // Stuff for calc
  double time_until_impact = p->getTimeLeft();
  double time_after_impact = timestep - p->getTimeLeft();

  glm::vec3 oldPos = p->getOldPos();
  glm::vec3 oldCen = p->getCenter();
  glm::vec3 dir = p->getDir();

  // We didn't hit an object in this interval of time
  if(time_until_impact >  timestep){

    glm::vec3 newPos( oldPos + dir * (float)timestep * VELOCITY_OF_MEDIUM );
    p->setPos(newPos);
    p->decTime(timestep);
  

  }else{
    // If we hit an object in this interval of time 

    // We where we hit in space
    glm::vec3 impactPos(oldPos + (dir * (float)time_until_impact) * VELOCITY_OF_MEDIUM);

    // Get the new center to change direction
    glm::vec3 mir_dir = MirrorDirection(p->getHitNorm(), p->getDir());
    mir_dir = mir_dir * (float)(-1.0);
    float radius = glm::distance(p->getCenter(), impactPos);

    // Move up to line
    p->setCenter(impactPos + mir_dir * radius);
    
    // p->setOldPos(impactPos);
    p->setPos(impactPos);

    calcMeshCollision(p); // new timeLeft Case a) we move a little up

    if(p->getMaterialHit() == "none"){

      // set center back to old center
      p->setCenter(impactPos + dir * radius);
      calcMeshCollision(p); // new timeLeft Case a) we move a little up
    
    }
 
    // Power has to be calculated after the hi:t
    double absorb_ratio = absorbFunc(p->getMaterialHit(), p->getFreq());
    assert( absorb_ratio < 1); // <----------------------------------------------------------------- Selfcheck
    p->setWatt( (1 - absorb_ratio ) * p->getWatt() ); // <------------------------------------------ Math 

    p->incIter();

    //moveParticle(p, time_after_impact); // this dude will move

    return true;

  }
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

 
    // Project back on sphere // When particles should die

    cirlce_point_on_sphere(p->getCenter(),
        glm::distance( p->getOldPos(), p->getCenter()),newPart);

    // For each calculated pos, make particle
    for(int i = 0; i < newPart.size(); i++){
    
      Particle * s = new Particle(
          newPart[i],                                   // Position
          newPart[i],                                   // OldPosition
          p->getCenter(),                               // CenterPos
          p->getWatt() / (double)(SPLIT_AMOUNT + 1.0),   // Amp
          p->getFreq(),                                 // Freq
          p->getSplit() + 1);                           // SplitAmount

      calcMeshCollision(s);                             // Manditory Calc
    
      // put particle there to be "moved" when its their turn
      vec.push_back(s);
    }// for
}

void ParticleSystem::calcMeshCollision(Particle * &p){
  // Input  :  Give a particle
  // Output :  None
  // Assumpt:  That direction and center are setup
  // SideEff:  Sets timeLeft until collision
  // SideEff:  Sets materialHit
  // SideEff:  Sets hitNormal
  
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

    if(triangle_intersect(r,h,a,b,c,backface)){
      hitTriangle = true;
      h.setMaterial(t->getMaterial());
    }
  }

  // Get the closest  hit
  if(hitTriangle){

    p->setTime(h.getT()/ VELOCITY_OF_MEDIUM);
    p->setHitNorm(h.getNormal());
    p->setMaterial(h.getMaterial());

  }else{

    // We don't hit any triangles, we most likey are a corner
    p->setTime(10000000); 
    p->setMaterial("none");

  }


}


void ParticleSystem::createInitWave(){
  // Testing function to create circle in 3d space

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

    // default constructor

    Particle * p = new Particle(
        pos,                     // Position     
        cursor,                  // OldPosition
        cursor,                  // CenterPos
        0,                       // Wattage
        0,                       // Freq
        0);                      // SplitAmount
  
    // For each given source type we will generate a diffrerent distrubution
    // of freq to represent
    
    if( args-> source_type == 1){

      // Low Freq Noise  AC unit
      // Not very powerfull

      p->setWatt(0.00001); // 70dBs about vacuum cleaner loudness
      p->setFreq(args->randomGen.randInt(80) + 20); // Freq [20Hz-100Hz]


    } else if (args-> source_type == 2 ){

      // Assume Talking Range
      // For speach we are using the fundmental, aka lowest freq
      
      p->setWatt(0.000001); // 60dBs about people talking loud
      p->setFreq(args->randomGen.randInt(150) + 100); // Freq [100Hz-250Hz]


    } else if( args-> source_type == 3 ){
      // Assume higher pitched noise
      // Assume you have a CRT mointor

      p->setWatt(0.0000001); // 40dDbs soft conversation level
      p->setFreq(16744); // Freq of CRT mointor running


    } else {

      // White noise there is total random distrubtion in frequency
      // Assume power of rock concert at 110 dBs

      p->setWatt(0.1); // Power of loud concert
      p->setFreq(args->randomGen.randInt(20000-20) + 20); // Freq random

    }
  
    // Find out when it hits our mesh
    calcMeshCollision(p);

    // put particle there
    particles.push_back(p);
  
  }
}

bool ParticleSystem::particleSplitCheckAndMerger(Particle * &p, std::vector<int> &deleteMask ){
  // Input:  A single particle, empty array filled by mask where 0 = keep 1 = merged and delete
  // Output: True if we should split on this particle, False otherwise.
  // Output: deleteMask will overwrite 1 at the index of particle that was merged
  // Assumptions: Particles concidered should exisit and not be merged
  // Assumptions: deleteMask the size of particle vector
  // Side Effects: Will replace the current particle with a new particle merges were required
  
  // FIXME: Using default split behavior 
  // ( split if there are not at least n particles some threhshold away from you )

  glm::vec3 pos = p->getOldPos();


  unsigned int particlesWithinTheshRequired = 3; // TODO should be global

  unsigned int particlesWithinThesh = 0;
  std::vector<Particle *> particleToMerge;

  float splitDistanceThresh = 5 * RADIUS_PARTICLE_WAVE;   // Play with these values
  float mergeDistanceThresh = 1 * RADIUS_PARTICLE_WAVE;
  float mergeAngleThesh = (2*M_PI) / 8.0;



  // for each particle in the system
  for(int i = 0; i < particles.size(); i++){

    // Do not count myself in this check and merger or particles already merged
    if(particles[i] == p || deleteMask[i] == 1)
      continue;

    // my distance to currently concidered particle
    float dist = glm::distance(p->getOldPos(), particles[i]->getOldPos());

    if(mergeDistanceThresh > dist ){
      // If we are close enough to concider merging these particles
      

      // Is the difference in direction similar
      float angle = acos( 
          glm::dot( particles[i]->getDir(), p->getDir() ) / 
          (glm::length(particles[i]->getDir()) * glm::length(p->getDir()))
      );
          
      if(mergeAngleThesh > angle){
      
        // Mark particle for deletion so we no longer concider it
        deleteMask[i] = 1;

        // Store for later merging
        particleToMerge.push_back(particles[i]);
      
      }
          
            

    }else if( splitDistanceThresh > dist) {
      // Still check if they fall within distance threshold for split check
      particlesWithinThesh++;
    
    }

  }//for


  // Handle Merges and overwrite this particle with a merged one
  if( particleToMerge.size() > 0){
  
    // Change that particle to new merged particle
    Particle * mergedPart =  particleVectorMerge(particleToMerge);
    *p = *mergedPart;
    return particleSplitCheckAndMerger(p, deleteMask);

  }
  
  // return false;
  return particlesWithinThesh < particlesWithinTheshRequired;

}


std::vector<std::vector<int>> ParticleSystem::createMunkresMatrix
  (const std::vector<Particle*> & partVec , const double angle){
  // Input: This function will take in a vector of particles to concider for
  //        the fitting a s hexagonal submask Also an angle for the 
  //        rotation of this submask
  //
  // Output: A square matrix where 1 represents a matching and 0 doesn't.
  //         
  //         
  //
  // Assumptions: The first element in the vector is the particle being
  //              concidered the center of the submask.
  //
  // Assumptions: There is at least a single element in the vector
  // Side Effects: None, works functionally
  //
  //
  // Output Example:
  //    
  //        MaskCenter MaskP1 MaskP2 MaskP3 MaskP4 MaskP5 MaskP6
  //  _______________________________________________
  //  partA   9         29     23      92     23    9
  //  partB   8         3      4       23      1    3 
  //  partC   9         1      3       1       90   3
  //  partD    .. .. ..
  //  .
  //  .
  //
  //
  //  Assumptions Example:
  //  [ a, b, c, d, e, f ...]
  //    ^
  //  This 'a' is concidered the center of my hexongal submask because it
  //  is first.
  //
  //
  //
  //  TODO: For reach rotation of the object, find the min cost and choose rotation
  //        of the best fit.
  //  NOTE: this version has no rotations applied


    // Checking we  have a valid input
    assert(partVec.size() != 0 );


    // ========================================================================
    // Create wrapper so that we can use the munkres.cpp
    // Remember to delete
    

    // Getting the center of our submask which we will base distance off of
    Particle * center = partVec[0];

    int sizeMatrix = 0; // 7 because 6 points + 1 center for hexagon
    partVec.size() > 7 ? sizeMatrix = partVec.size() : sizeMatrix = 7;

    // Create a square matrix and initiate all of it with infinate values
    int ** matrix = new int*[sizeMatrix];

    for(int i = 0; i < sizeMatrix; i++){

       // create a new row
       matrix[i] = new int[sizeMatrix];

       // fill in with super large value
       for( int j = 0; j < sizeMatrix; j++){
         matrix[i][j] = std::numeric_limits<int>::max();
       
       }
    }
    

    // Getting the points that represent the hyprotheical mask
    std::vector<glm::vec3> maskPositions; 
    circle_points_on_plane(center->getOldPos(), center->getDir(), RADIUS_PARTICLE_WAVE,6,maskPositions,args);


    // Comparing the distance between each point in partVec and each point in
    // point in the submask and tossing it into the matrix. 
    // Note: there will be partVec.size() * (6 + 1) comparisons

    // For every point
    for(int i = 0; i < partVec.size(); i++){

      // For every possible mask postion
      for(int j = 0; j < maskPositions.size(); j++){

        glm::vec3 partPos = partVec[i]->getOldPos();

        // distance calculation to get cost
        double dist = glm::distance(partPos, maskPositions[j]);

        // This will put in the scale of milimeters everything inside my matrix
        // Of which is small enough scale that it cover high freq wave lengths
        matrix[i][j] = (int) (dist * 1000);
      
      }
    }

  // Call to munkers.cpp hungarian algorithm
  matrix = runMunkers(matrix,sizeMatrix,false);

  // Conversion  back into vector
  std::vector<std::vector<int>> matching;
  for(int i = 0; i < partVec.size(); i++){
    
    std::vector<int> temp;
    for(int j = 0; j < maskPositions.size(); j++){
      temp.push_back(matrix[i][j]);
    }
    matching.push_back(temp);
  }
    
  // deallocation of created array in matrix **
  for(int i = 0; i < sizeMatrix; i++)
    delete [] matrix[i];
  delete [] matrix;
  
  // Return the matching for use
  return matching;

}


bool ParticleSystem::shouldSplit(Particle * &p){
  // In here we compare all particles against eachother
  // Very expensive to do, however we also check to see if
  // we can merge any two particles into one
  return false; //<-------------------------------------------------------------DEBUG

  glm::vec3 pos = p->getOldPos();
  float nearestDistance = 100000;
  float threshold = 3 * RADIUS_PARTICLE_WAVE;

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

Particle * ParticleSystem::particlePairMerge(Particle * &a, Particle * &b){ 

  // Input: A pair of particles and an empty pointer that will be result in a new particle
  // Output: None, handled through pass by reference
  // Assumptions: non null pointers
  // Side Effects: Calls directly particleVectorMerge


  std::vector<Particle *> vec;
  vec.push_back(a); vec.push_back(b);

  return particleVectorMerge(vec);

}


Particle * ParticleSystem::particleVectorMerge(std::vector<Particle *> &vec){

  // Input: A vector of particles and an empty pointer that will be result in a new particle
  // Output: None, handled through pass by reference
  // Assumptions: vec is populated with at least 1 particle 
  // Side Effects: None
  // Bugs: What should I do with the freq of merged particles & splits (split not as much)

  // Values will accumlate to find the average
  double pos_x = 0; double cen_x = 0;
  double pos_y = 0; double cen_y = 0;
  double pos_z = 0; double cen_z = 0;


  double mergeFreq = 0; double mergeWatt = 0;
  double mergeSplit = 0;
  int    mergeIter = 0;

  for( int i = 0;  i <  vec.size();  i++ ){

    const glm::vec3 oldPos = vec[i]->getOldPos();
    const glm::vec3 cen = vec[i]->getCenter();

    pos_x += oldPos.x; cen_x += cen.x; 
    pos_y += oldPos.y; cen_y += cen.y;
    pos_z += oldPos.z; cen_z += cen.z;


    mergeFreq  += vec[i]->getFreq();
    mergeWatt  += vec[i]->getWatt();
    mergeSplit += vec[i]->getSplit();


    // Find max iterations among particles in vec
    if( vec[i]->getIter() > mergeIter )
      mergeIter = vec[i]->getIter();

  }

  glm::vec3 mergedPos(
      pos_x / vec.size(), 
      pos_y / vec.size(), 
      pos_z / vec.size());

  glm::vec3 mergedCen(
      cen_x / vec.size(), 
      cen_y / vec.size(), 
      cen_z / vec.size());

  mergeFreq  = mergeFreq  / vec.size(); //FIXME What should be done regarding freq ?
  mergeSplit = mergeSplit / vec.size();

  Particle * newPart = new Particle(mergedPos,mergedPos,mergedCen, mergeWatt, mergeFreq, mergeSplit);
  newPart->setIter(mergeIter);

  // Calc where I will hit next
  calcMeshCollision(newPart);


  return newPart;

}

double ParticleSystem::absorbFunc(const std::string & mtlName, 
    const double freq){

  std::string materialName = mtlName; // <---- change back to materialName 

  // A few checks
  if( freq < 20 ){
    std::cout << "Freq too low" << std::endl;
    assert(false);
  }


  // Remapper 
  
  if( materialName.compare(0,5,"GLASS") == 0){
    materialName = "double_window";
  }

  // We will assume only some walls is made of Parete Absorber
  else if(materialName.compare(0,4,"wall") == 0){

      // Get the index
      std::string temp  = materialName.substr(5);
      int index = atoi(temp.c_str()) % 3;

      if(index == 1 && args->absorber) { //<------------------------------------------ Get absorber to show

        // Make this wall an absorber
        materialName = "absorber_parete";

      }else{

        // Make these walls bricked walls
        switch (args->wall_material) {
          case 0:
            materialName = "bricked_wall";
            break;
          case 1:
            materialName = "concrete_wall";
            break;
          case 2:
            materialName =  "ceramic_wall";
            break;
          default:
            assert(false);
        }

      }
  } // walls

  else if( materialName.compare(0,6,"FILLIN") == 0) {
  
    switch (args->wall_material) {
      case 0:
        materialName = "bricked_wall";
        break;
      case 1:
        materialName = "concrete_wall";
        break;
      case 2:
        materialName =  "ceramic_wall";
        break;
      default:
        assert(false);
    }
  
  }

  else if( materialName == "floor"){
  
    switch (args->floor_material) {
      case 0:
        materialName = "pvc_floor";
        break;
      case 1:
        materialName =  "carpated_floor";
        break;
      default:
        assert(false);
    }
  
  }

  else{
  
    // Assume ceiling
    materialName = "plaster_ceiling";
  
  }

  // Name remapper ////////////////////////////////////////////////////////////
  
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
      return 0.00002625 * (freq - 507) + 0.01;

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

