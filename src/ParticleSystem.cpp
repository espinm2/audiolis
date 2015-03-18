#include "ParticleSystem.h"

#include <glm/glm.hpp>
#include <glm/gtx/vector_angle.hpp>
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
#include "mask.h"
#include <algorithm>


#define EPSILON 0.0001
typedef unsigned int uint;
typedef short unsigned int uint8;
typedef std::vector<Particle *>  PartPtrVec;

// Used for update
typedef std::vector<Particle *>::iterator ParticleIter;
typedef std::vector<std::vector<int>> vMat;

ParticleSystem::~ParticleSystem(){
  // Just delete all the particles we made
  for(unsigned int i = 0; i < particles.size(); i++)
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
  VELOCITY_OF_MEDIUM    = 340; //implement to sound speed in m/s
  if(args->timestep != -1 )
    TIME_STEP = args->timestep ;
  else
    TIME_STEP = 0.001 / VELOCITY_OF_MEDIUM; // Centemimer accruacy

  // simuation var init
  RADIUS_INIT_SPHERE    = 0.7; // Doesnt get changed to ofte
  NUM_INIT_PARTICLES    = args->num_init_particles; // imported from args
  MIN_WATTAGE           = 0.000000000002; // Doesnt change often
  MAX_ITERATIONS        = 6000; // Doesnt change often

  // split var init
  RADIUS_PARTICLE_WAVE  = 0.1; // Will be later removed for better split
  SPLIT_AMOUNT          = 6; // Will be later removed for better split

  ITERATION             = 0;

  // Profiler stuff
  if(args->profile){
    output_profiler_str.open(args->output_file);
    if( !output_profiler_str.good() ){
      std::cerr << "Cant open output file for profiler " << std::endl;
      exit(1); }
 }
}

void ParticleSystem::update(){

  // Clear the new positions
  for(Particle * cur : particles){
    cur->setOldPos(cur->getPos());
    cur->setPos((glm::vec3)NULL);
  }

  // Move all the partilces in the system up a timestep
  for(Particle * cur : particles)
    moveParticle(cur);

  // // Create all possible splits
  // for(Particle * cur : particles)
  //   generateResSplits(cur);

  // // New particles are added to the mix
  // for(Particle * np : newParticles)
  //   particles.push_back(np);

  // // Merge particles that are the same
  // for(Particle * cur : particles)
  //   mergeSimilarParticles(cur);

  // Resolve all collision that occur
  for(Particle * cur: particles)
    resolveCollisions(cur);

  // // Cleans our particles vector and our newParticles
  // removeDeadParticles();
  // newParticles.clear();

  // Remakes the kd tree
  // particle_kdtree.update(particles, *bbox);
}

void ParticleSystem::moveParticle(Particle* &p){
 /* Input : Particle ptr
  * Output: None
  * Asumpt: That particle exisit
  * SideEf: Changes position of the particle & up incr
  */

  assert(p != NULL && p->isAlive());

  // Calculate new position
  glm::vec3 oldPos = p->getOldPos();
  glm::vec3 dir = p->getDir();
  glm::vec3 newPos( oldPos + dir * (float)TIME_STEP * VELOCITY_OF_MEDIUM );

  // Update the particle
  p->setPos(newPos); p->incIter();
}

void ParticleSystem::resolveCollisions(Particle* &p){

  // Create a ray & hit class
  Ray r(p->getOldPos(), p->getDir());
  Hit h; bool hitTriangle; bool backface = false;

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

  // If we hit nothing, we are out in the open
  if(!hitTriangle){ std::cout << "PS: ERROR_OUT_OF_BOUNDS\n"; return; } 

  // If we hit do not hit within our timestep
  if(h.getT() + EPSILON > TIME_STEP){ return; }

  // Gareneteee that we hit just this timestep
  // Change the direction of our particle & backstep to hitting wall
  double time_until_impact = h.getT();
  glm::vec3 old = p->getOldPos();
  glm::vec3 dir = p->getDir();
  glm::vec3 impactPos(old+(dir*(float)time_until_impact)*VELOCITY_OF_MEDIUM);

  // Get the new center to change direction
  glm::vec3 mir_dir = MirrorDirection( h.getNormal() , p->getDir() );
  mir_dir = mir_dir * (float)(-1.0);
  float radius = glm::distance(p->getCenter(), impactPos);

  // Mirror the center to reflect our new direction
  p->setCenter(impactPos + mir_dir * radius);

  // stick directly into the wall (max loosing 1mm)
  p->setPos(impactPos);


}

void ParticleSystem::generateResSplits(Particle * &cur){
 /* Input : Particle ptr
  * Output: None
  * Asumpt: That particle exisit
  * SideEf: Adds to new particle vector
  */
  assert(cur!= NULL && cur->isAlive());

  // Properties we will use for gathering and merging
  double gather_distance = RADIUS_PARTICLE_WAVE * 2.5;
  double gather_angle    = M_PI / 16.0; 

  // temprary vectors
  PartPtrVec gathered_particles;           // temp particle holders
  std::vector < glm::vec3 > new_positions; // used incase we split

  // requirment for mask fitting
  gathered_particles.push_back(cur);

  // Gather the particles I will use for mask fitting
  particle_kdtree.GatherParticles(cur, 
  gather_distance, gather_angle, gathered_particles);

  // Generate my mask
  Mask mask;
  generateMask(gathered_particles, mask);

  if( !mask.resSpit( new_positions )){ return; } // skip no splits occur

  // Create new particles
  for(glm::vec3 pos : new_positions){

    Particle * s = new Particle(
      pos,  // WARNING THIS MIGHT NEED TO CHANGE
      pos, 
      cur->getCenter(), 
      cur->getWatt() / (double)(SPLIT_AMOUNT + 1.0), 
      cur->getFreq(), 
      cur->getSplit() + 1 
    );

    newParticles.push_back(s);

  }
}

void ParticleSystem::mergeSimilarParticles(Particle * &cur){
  // Implement
  assert(false);
}

void ParticleSystem::stabalizeInitalSphere(){
  // TODO: Test this function

  // Where we store how much force is in the inital system
  glm::vec3 sumForces(0,0,0);

  for(Particle * p : particles){

    p->setOldPos(p->getPos());
    // Find nearest particle
    float nearest_distance = 1000;
    glm::vec3 nearest_pos = (glm::vec3) NULL;

    for(Particle * o : particles) {

      if( o == p ) // dont count self
        continue;

      float dist = glm::distance(p->getOldPos(), o->getOldPos());

      // Are we close enough
      if( dist < nearest_distance ){
        nearest_distance = dist;
        nearest_pos = o->getOldPos();
      }
    }

    // Case where no particles are near me 
    if(nearest_distance == 1000)
      continue;
    
    // Get force
    glm::vec3 force = CalcRepulsiveForces(
      p->getOldPos(), nearest_pos, RADIUS_PARTICLE_WAVE, 0.05);

    sumForces += force;

    // Apply force to my particle x0 + 0.5 f t^2
    glm::vec3 newPos = p->getOldPos()  +  force;

    newPos = glm::normalize(newPos - p->getCenter()) 
    * (float)RADIUS_INIT_SPHERE + p->getCenter();
    
    p->setPos(newPos);

  }

  // Check if we are stabalized enough yet
  if(glm::length(sumForces) < 0.007){
    args->setupInitParticles = false;
    std::cout << "finished stabalizeInitalSphere: " << glm::length(sumForces) << "\n";
  }
}

void ParticleSystem::linearGatherParticles(Particle * center, double r, double a, PartPtrVec & result){

  // Old Gather code
  for(unsigned int i = 0; i < particles.size(); i++){

    Particle * other = particles[i];

    if(other == center)
      continue;

    if(other->isDead())
      continue;

    float dist = glm::distance(center->getOldPos(), other->getOldPos());

    if(dist < r){

      float angle = angleBetweenVectors(center->getDir(), other->getDir());

      if( angle < a)
        result.push_back(other);
    }
  }
}

bool ParticleSystem::linearDuplicateSearch(const glm::vec3 & pos, double th){
   for( Particle * p : particles)
     if ( glm::distance(p->getOldPos(), pos ) < th)
      return true;
  return false;
}

bool ParticleSystem::linearNewDuplicateSearch(const glm::vec3 & pos, const PartPtrVec & newVec , double th){

   for( Particle * p : newVec)
     if ( glm::distance(p->getOldPos(), pos ) < th)
      return true;
  return false;
}


void ParticleSystem::moveCursor( const float & dx, 
    const float & dy, const float & dz ){
  // Function called by glCanvas
  
  cursor+= glm::vec3(10*dx,10*dy,10*dz);
  // std::cout << cursor << std::endl;
}

void ParticleSystem::createDebugParticle(){

  // HARDCODED targetPosition
  glm::vec3 targetPosition(2.7, 0.324, -2.1); // for corner room
  // glm::vec3 targetPosition(-5.8, 1.524, -0.2); // for acoustics 

  // Direction 
  glm::vec3 directionToTarget = targetPosition - cursor;
  directionToTarget  = glm::normalize(directionToTarget);

  // Create a ray to move up a little
  glm::vec3 newPos = 
    cursor + ( (float) RADIUS_INIT_SPHERE * 2  ) * directionToTarget;



  // default constructor
  Particle * p = new Particle(
      newPos,                     // Position     
      newPos,                  // OldPosition
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

  // Get new particles to be made and toss in split particles
  std::vector<Particle *> new_particles;
  particleSplit(p, new_particles);
  new_particles.push_back(p);

  // change cur particle watts
  p->setWatt(p->getWatt() / (double)(SPLIT_AMOUNT + 1.0));
  
  for(Particle * sp: new_particles)
    particles.push_back(sp);

  particle_kdtree.update(particles,*bbox);

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

 
    // Project back on sphere // When particles should die

    // cirlce_point_on_sphere(p->getCenter(),
    //    glm::distance( p->getOldPos(), p->getCenter()),newPart);

    // For each calculated pos, make particle
    for(unsigned int i = 0; i < newPart.size(); i++){
    
      Particle * s = new Particle(
          newPart[i],                                   // Position
          newPart[i],                                   // OldPosition
          p->getCenter(),                               // CenterPos
          p->getWatt() / (double)(SPLIT_AMOUNT + 1.0),   // Amp
          p->getFreq(),                                 // Freq
          p->getSplit() + 1);                           // SplitAmount

    
      // put particle there to be "moved" when its their turn
      vec.push_back(s);
    }// for
}


void ParticleSystem::createInitWave(){
  // Testing function to create circle in 3d space

  // Using cursor
  double s = RADIUS_INIT_SPHERE;
  
  // Create Box of ranodm points
  for(unsigned int i = 0; i < NUM_INIT_PARTICLES; i++){
    // Find x,y,z
    float x = cursor.x - s/2.0 + (float) args->randomGen.rand(s);
    float y = cursor.y - s/2.0 + (float) args->randomGen.rand(s);
    float z = cursor.z - s/2.0 + (float) args->randomGen.rand(s);
  
    glm::vec3 pos(x,y,z);
  
    // Project into a circle
    float radius = s;
    
    glm::vec3 dir = pos - cursor;
  
    dir = glm::normalize(dir);
  
    pos = cursor + dir * radius;

    // default constructor

    Particle * p = new Particle(
        pos,                     // Position     
        pos,                     // OldPosition
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
  
    // put particle there
    particles.push_back(p);
  
  }


  particle_kdtree.update(particles, *bbox);
}


void ParticleSystem::munkresMatching 
  (std::vector<Particle*> & partVec, vMat & matchingMat, vMat & costMat){

  // Function: Creates a matching to a hexonagal mask on particles given
  //
  // Input: partVec is a function of particles, partVec[0] will be the center
  //        of this mask.
  //
  // Output: A 2D vector where 1 == matching and 0 == not matched.
  //         Column = specific mask postion we are matching against
  //         Row    = Row i, means particle i in our partVec
  //         
  // Assumptions: The first element in the vector is the particle being
  //              concidered the center of the submask.
  //
  // Assumptions: There is at least a single element in the vector
  //
  // Side Effects: None.
  //
  //
  // Output Example:
  //    
  //        MaskCenter MaskP1 MaskP2 MaskP3 MaskP4 MaskP5 MaskP6
  //  _______________________________________________
  //  partA   1         0      0       0       0     0      0
  //  partB   0         1      0       0       0     0      0
  //  partC   0         0      0       1       0     0      0
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
  

  assert(partVec.size() != 0 ); // Checking we  have a valid input 

  // ==========================================================================
  // Create wrapper matrix so that we can use the munkres.cpp
  // ==========================================================================
  
  // Getting the center of our submask which all distance will be relative to
  Particle * center = partVec[0];
  glm::vec3 centerMaskPos = center->getOldPos();

  // Setting the dimension of our square matrix
  unsigned int maxtrix_dimension = 0;
  partVec.size() > 7 ? maxtrix_dimension = partVec.size() : maxtrix_dimension = 7;

  // Create a square matrix and initiate all of it with infinate values
  int ** matrix = new int*[maxtrix_dimension];
  for(unsigned int i = 0; i < maxtrix_dimension; i++){
     // create a new row
     matrix[i] = new int[maxtrix_dimension];
     // fill in with super large value
     for(unsigned int j = 0; j < maxtrix_dimension; j++){
       matrix[i][j] = 1000; // This is an entire 1000mm = 1m 
     }
  }
  
  // ==========================================================================
  // Generation of delusional particle positions
  // ==========================================================================

  // Get the nearest particle around the center of the mask
  // We will use this particle to orient our delusional particles
  
  // Where we store the points that represent the hyprotheical mask
  std::vector<glm::vec3> maskPositions; 
  maskPositions.push_back(centerMaskPos); // maskPositions[0] is mask center


  delusionalParticleLocations(center,partVec, maskPositions);

  // ==========================================================================
  // Creating a cost for each particle being concidered
  // ==========================================================================
  
  // Comparing the distance between each point in partVec and each point in
  // point in the mask and tossing it into the matrix. 
  // Note: there will be partVec.size() * (6 + 1) comparisons

  // For every point
  for(unsigned int i = 0; i < partVec.size(); i++){

    // Concidering my cur, we should only let it match with itself
    if(i == 0){
      matrix[0][0] = 0;
      continue;
    }

    // For every possible mask postion
    for(unsigned int j = 0; j < maskPositions.size(); j++){

      // Position of particle I am concidering right now
      glm::vec3 partPos = partVec[i]->getOldPos();

      // distance calculation from partPos to all possible mask
      double dist_from_ideal = glm::distance(partPos, maskPositions[j]);
      bool within_angle_constrant = true;

      glm::vec3 dir_ideal = 
        glm::normalize(maskPositions[j] - center->getOldPos());

      glm::vec3 dir_of_partPos = 
        glm::normalize(partPos-center->getOldPos());

      double direction_dist = glm::distance(dir_of_partPos,dir_ideal);

      if( dist_from_ideal > EPSILON && direction_dist > EPSILON && j != 0){ // conditions that break angle calc

        double angle = acos( glm::dot( dir_ideal, dir_of_partPos ) / 
          (glm::length( dir_ideal ) * glm::length( dir_of_partPos)));

        angle <= 0.436332313 ? within_angle_constrant = true : within_angle_constrant = false;

      }

      // Prevent concave shapes
      if( dist_from_ideal <= 1.2*RADIUS_PARTICLE_WAVE && within_angle_constrant  ){ 
    

        // This will put in the scale of milimeters everything inside my matrix
        matrix[i][j] = (int) (dist_from_ideal * 1000); 
    
      }
    }
  }


  // Conversion from 2d array to 2d vector
  for(unsigned int i = 0; i < partVec.size(); i++){

    std::vector<int> temp;

    for(unsigned int j = 0; j < maskPositions.size(); j++)

      temp.push_back(matrix[i][j]);

    costMat.push_back(temp);

  }


  // ==========================================================================
  // Compute the optimal solution
  // ==========================================================================
  matrix = runMunkers(matrix,maxtrix_dimension,false);


  // ==========================================================================
  // Save the result
  // ==========================================================================
  for(unsigned int i = 0; i < partVec.size(); i++){

    std::vector<int> temp;

    for(unsigned int j = 0; j < maskPositions.size(); j++)
      temp.push_back(matrix[i][j]);

    matchingMat.push_back(temp);
  }

  
  // ==========================================================================
  // deallocation of created array in matrix **
  // ==========================================================================
  
  for(unsigned int i = 0; i < maxtrix_dimension; i++)
    delete [] matrix[i];
  delete [] matrix;
}


void ParticleSystem::generateMask(
    std::vector <Particle*> & conciderForMask, Mask &m ){
  // Input: p is a particle that it not null, it will be the center of the mask
  // Input: Mask &m is a mask that will be populated by a mask object
  // Output: None, see input (2)
  
  // Assumptions: conciderForMask[0] is the center of the mask
  // Side effects: None.
  // Preformance: O(p) where p is the size of particles vector


  vMat cost; vMat matching;

  // Calculate the cost and matching matries
  munkresMatching(conciderForMask,matching,cost);


  /*
  // Uncomment for debugging cost matrix 
  std::cout << "Cost Matrix \n";

  int p_index = 0;
  for(std::vector<int> v : cost){
    std::cout << &(*conciderForMask[p_index]) << "\t[ ";
    for(int c : v){
      std::cout << c << "\t";

    }
    std::cout << "]\n";
    p_index++;
  }

  std::cout << "Matching Matrix \n";

  p_index = 0;
  for(std::vector<int> v : matching){
    std::cout << &(*conciderForMask[p_index]) << "\t[ ";
    for(int c : v){
      std::cout << c << "\t";

    }
    std::cout << "]\n";
    p_index++;
  }
  // Uncomment for debugging cost matrix end */

  assert(matching.size() == cost.size());
  assert(matching[0].size() == cost[0].size());

  assert(matching.size() == conciderForMask.size()); 
  assert(matching[0].size() == 7); 

  // Inners of the mask class
  std::vector<Particle*> maskPart;
  std::vector<int> maskCost;
  int size_of_mask = 0;


  // For each column in the matching matrix push back the matching particle 
  for(unsigned int j = 1; j < matching[0].size(); j++){
    
    bool found = false;

    // Go and find what particle matches this
    for(unsigned int i = 1; i < matching.size(); i++){

      if(matching[i][j] == 1 && cost[i][j] < 1000){ // second part is because this makes arb matches

        maskPart.push_back(conciderForMask[i]);
        maskCost.push_back(cost[i][j]);
        size_of_mask++;
        found = true;
        break; // found particle from column
      }
    }

    if(!found){
      // If we can't find a match push back null
      maskPart.push_back(NULL);
      maskCost.push_back(1000); // cost is 1000mm = 1m
    }

  }// for each column

  assert(maskPart.size()==6);
  // Setting the memebers of the mask object
  m.setCenter(conciderForMask[0]);
  m.setMaskParticles(maskPart);
  m.setCostVector(maskCost);
  m.setSize(size_of_mask);
}

Particle * ParticleSystem::particlePairMerge(Particle * &a, Particle * &b){ 

  // Input: A pair of particles and an empty pointer that 
  //   will be result in a new particle
  // Output: None, handled through pass by reference
  // Assumptions: non null pointers
  // Side Effects: Calls directly particleVectorMerge


  std::vector<Particle *> vec;
  vec.push_back(a); vec.push_back(b);

  return particleVectorMerge(vec);
}


Particle * ParticleSystem::particleVectorMerge(std::vector<Particle *> &vec){

  // Input: A vector of particles and an empty pointer that will be 
  //  result in a new particle
  // Output: None, handled through pass by reference
  // Assumptions: vec is populated with at least 1 particle 
  // Side Effects: None
  // Bugs: What should I do with the freq of merged particles & splits 
  //    (split not as much)

  // Values will accumlate to find the average
  double pos_x = 0; double cen_x = 0;
  double pos_y = 0; double cen_y = 0;
  double pos_z = 0; double cen_z = 0;


  double mergeFreq = 0; double mergeWatt = 0;
  double mergeSplit = 0;
  int    mergeIter = 0;

  for(unsigned int i = 0;  i <  vec.size();  i++ ){

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

  mergeFreq  = mergeFreq  / vec.size();                                         // Patch for merging freq
  mergeSplit = mergeSplit / vec.size();

  Particle * newPart = new Particle(
      mergedPos,mergedPos,mergedCen, mergeWatt, mergeFreq, mergeSplit);
  newPart->setIter(mergeIter);



  return newPart;
}

double ParticleSystem::absorbFunc(const std::string & mtlName, 
    const double freq){

  std::string materialName = mtlName;

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

      if(index == 1 && args->absorber) { 

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

void ParticleSystem::delusionalParticleLocations(
    Particle * &cur_particle,
    std::vector<Particle *> &gathered_particles,
    std::vector<glm::vec3> & output){

  // Input : cur_particle is the particle we will use the center of the mask
  // Input : gathered_particles are the particles we will use to orient ourselves
  // Input (output) :  We will return our locations here


  assert(cur_particle == gathered_particles[0]);


  // Search for neareast particle
  double dist = RADIUS_PARTICLE_WAVE * 10; // really large number
  unsigned int nearest_part = 0;

  for(unsigned int i = 1; i < gathered_particles.size(); i++){
    float  curDist = glm::distance(gathered_particles[i]->getOldPos(), cur_particle->getOldPos());
    if( curDist < dist  ){
        dist = curDist; nearest_part = i;
      }
  }

  glm::vec3 nearest_pos = gathered_particles[nearest_part]->getOldPos();
  
  // We will adjust the mask size of our particles by a factor of two
  // If we are within the range of [0, 2*RADIUS_PARTICLE_WAVE]
  double mask_radius = dist;

  // Cap our distance at 
  if( mask_radius > 2.2 * RADIUS_PARTICLE_WAVE )
    mask_radius = 2.2 * RADIUS_PARTICLE_WAVE;


  circle_points_on_plane_refence(
      cur_particle->getOldPos(),                        // center of mask
      cur_particle->getDir(),                     // direction of plane
      nearest_pos,              // particle using for reference
      mask_radius,                 // radius of my mask
      6,                                    // number of particles mask has
      output);                       // where I will append my results


  // We project there mask particles  on the sphere of the sound source
  circle_points_on_sphere(
    cur_particle->getCenter(),                                      // center
    glm::distance(cur_particle->getCenter(), cur_particle->getOldPos()), // radi
    output);


  /*
circle_points_on_sphere
  assert(cur_particle == gathered_particles[0]);
  circle_points_on_plane(
      cur_particle->getOldPos(),                        // center of mask
      cur_particle->getDir(),                     // direction of plane
      // nearestParticlePosition,              // particle using for reference
      RADIUS_PARTICLE_WAVE,                 // radius of my mask
      6,                                    // number of particles mask has
      output);                       // where I will append my results


  */
}
