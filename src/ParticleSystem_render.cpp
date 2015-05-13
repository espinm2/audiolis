#include <glm/glm.hpp>
#include <math.h>

#include "ParticleSystem.h"
#include "particle.h"
#include "boundingbox.h"
#include "camera.h"
#include "argparser.h"
#include "render_utils.h"
#include "geometry_utils.h"
#include "sphere.h"

#define MAX_ITERATIONS 6000
#define MAX_ITERATIONS 6000
#define MIN_WATTAGE 0.000000000002

// This is for rainbow color
#define BLUE     getColor(1,66,255,1)
#define LITE_BLU getColor(0,254,252,1)
#define GREEN    getColor(135,255,119,1)
#define YELLOW   getColor(255,240,1,1)
#define RED      getColor(255,30,0,1)

#define OUTLINE_VIZ true

typedef std::vector<std::vector<int>> vMat;

void ParticleSystem::initializeVBOs(){

  // Get uniquie ids from buffers
  glGenBuffers(1,&particle_verts_VBO);
  glGenBuffers(1,&cursor_verts_VBO);
  glGenBuffers(1,&velocity_verts_VBO);
  glGenBuffers(1,&velocity_tri_indices_VBO);
  glGenBuffers(1,&outline_verts_VBO);
  glGenBuffers(1,&happyness_verts_VBO);
  glGenBuffers(1,&delusional_verts_VBO);
  glGenBuffers(1,&connection_verts_VBO);
  glGenBuffers(1,&sphere_verts_VBO);
  glGenBuffers(1,&sphere_tri_indices_VBO);

  particle_kdtree.initializeVBOs();
  uniform_grid.initializeVBOs(); 
}

void ParticleSystem::setupVBOs(){
  HandleGLError("enter setup vbos");

  // Delete old data
  cursor_verts.clear();
  particle_verts.clear();
  velocity_verts.clear();
  velocity_tri_indices.clear();
  outline_verts.clear();
  happyness_verts.clear();
  delusional_verts.clear();
  connection_verts.clear();
  sphere_verts.clear();            // verts
  sphere_tri_indices.clear();
  
  // Setup new Data
  if(args->kdtree_render)
    particle_kdtree.setupVBOs();

  if(args->ugrid_render)
    uniform_grid.setupVBOs();

  setupCursorPoint();

  if(args->direction)
    setupVelocityVisual();

  if(args->viz_type==4)
    setupSphere();

  setupParticles(); 

  if(args->render_edges){
    // setupEdges(); // Dead to me
    setupDelusionalParticles(); // better visualization
  }

  HandleGLError("leave setup vbos");
}

void ParticleSystem::drawVBOs(){
  HandleGLError("enter draw vbos");

  drawCursorPoint();

  if(args->render_edges){
    // drawHappinessVisual();
    drawDelusionalParticles();
    drawDelusionalConnections();
  }

  if(args->direction)
    drawVelocityVisual();

  if(args->kdtree_render)
    particle_kdtree.drawVBOs();

  if(args->ugrid_render)
    uniform_grid.drawVBOs();

  if(args->viz_type == 4)
    drawSphere();

  drawParticles();

  HandleGLError("leaving draw vbos");

}

void ParticleSystem::cleanupVBOs(){
  HandleGLError("enter cleanup vbos");
  glDeleteBuffers(1,&particle_verts_VBO);
  glDeleteBuffers(1,&cursor_verts_VBO);
  glDeleteBuffers(1,&velocity_verts_VBO);
  glDeleteBuffers(1,&velocity_tri_indices_VBO);
  glDeleteBuffers(1,&outline_verts_VBO);
  glDeleteBuffers(1,&happyness_verts_VBO);
  glDeleteBuffers(1,&delusional_verts_VBO);
  glDeleteBuffers(1,&connection_verts_VBO);
  glDeleteBuffers(1,&sphere_verts_VBO);
  glDeleteBuffers(1,&sphere_tri_indices_VBO);

  HandleGLError("Space clean enter");
  particle_kdtree.cleanupVBOs();
  uniform_grid.cleanupVBOs();
  HandleGLError("space clean leave");


  HandleGLError("leave clean vbos");
}


uint ParticleSystem::getLowestCostShape(Particle * cur){
  // inputs:
  //      cur is the particle we are looking to find the best fit around it
  // output:
  //      uint a number from 0,[4:7]

  // Gather the particles I will use to make an analysis of points 
  PartPtrVec gathered_particles;        
  particle_kdtree.GatherParticles(
    cur,                              // This particle we collect around
    RADIUS_PARTICLE_WAVE * 3.0,       // collect a bit beyond mask
    GATHER_ANGLE,                     // if angle beyond this do not gather
    gathered_particles);              // where we will place these particles


  // Not enough particles to run, just return 0
  if(gathered_particles.size() < 4 )
    return 0;

  // Populate scoreVector
  std::vector <double > scoreVector;
  for( uint sides_shape = 4; sides_shape < 10; sides_shape++){

    if(gathered_particles.size() < sides_shape ){
      scoreVector.push_back(std::numeric_limits<double>::max());
      continue;
    }

    // Find the nearest 6 (sort by distance to center particle)
    const glm::vec3 centerPos = cur->getOldPos();

    std::sort(
      gathered_particles.begin(), 
      gathered_particles.end(), 
      particleCMP(centerPos)
    );

    PartPtrVec nearbyInOrder(
      gathered_particles.begin(), 
      gathered_particles.begin() + sides_shape
    );

    assert(nearbyInOrder.size() == sides_shape );

    // Pick a refrence particle
    Particle * ref = nearbyInOrder[0];

    std::vector<double> angles;
    for( uint i = 0; i < nearbyInOrder.size(); i++ ){

      double absAngle =                                   
      getAbsAngle(
        ref->getOldPos(),
        nearbyInOrder[i]->getOldPos(),
        cur->getOldPos(),
        glm::normalize(cur->getOldPos() - cur->getCenter())
      );

      angles.push_back(absAngle);

    }

    // Sort by angle 
    std::sort(angles.begin(),angles.end());

    // Get the squared error
    std::vector<double> error; // where i will keep all my error
    for(uint i = 0; i < angles.size(); i++){
      double deg = (360 / sides_shape) * i;
      error.push_back(pow( deg - angles[i]  ,2));
    }

    double error_total = 0;
    for(double e: error)
      error_total += e;
    error_total /= sides_shape;

    scoreVector.push_back(error_total);
  }// for shape size n

  unsigned int best_shape_fit = -1;
  double min_score = std::numeric_limits<double>::max();

  // Finding Max
  for(int i = 0; i < scoreVector.size(); i++){

    if(min_score > scoreVector[i]){
      best_shape_fit = i;
      min_score = scoreVector[i];
    }
  }

  assert(best_shape_fit != -1);

  return best_shape_fit + 4;

}

void ParticleSystem::setupParticles(){

  // Push everything into particle_verts
  for(int i = 0 ; i < particles.size(); i++){

    // Cur Particle
    Particle * part = particles[i];

    // Getting pos (current)
    glm::vec3 pos = part->getPos();

    // Picking Normal (normal direction)
    glm::vec3 normal =part->getDir(); 

    
    glm::vec4 color(0,0,0,1);  // Our final color


    if( args->viz_type == 0 ) {

      glm::vec4 colorA(0,0,0,1); // White -- Low Freq
      glm::vec4 colorB(1,1,1,1); // Black -- High Freq
    
      // Frequency visualization
      double freq = part->getFreq();
      
      assert( 20 <= freq && freq <= 20000); // Error if out of range
      double val = ( freq - 20 ) / ( 20000 - 20 ); // part/whole

      color.r = colorA.r + val * (colorB.r - colorA.r);
      color.g = colorA.g + val * (colorB.g - colorA.g);
      color.b = colorA.b + val * (colorB.b - colorA.b);
    } 
    
    else if (args->viz_type == 1){

      glm::vec4 colorA(1,1,1,1); // Wihte -- Low inseinty (runnig)
      glm::vec4 colorB(0,0,0,1); // Black -- lots of energy

      // Convert into dBs
      double dbs = 10 * log10( part->getWatt() / MIN_WATTAGE);
      double val = dbs / 120.0;


      // Sanity check
      if(dbs < 0 || dbs > 120){
        //std::cout << dbs << std::endl;
        val = 0.0;
      }

      color.r = colorA.r + val * (colorB.r - colorA.r);
      color.g = colorA.g + val * (colorB.g - colorA.g);
      color.b = colorA.b + val * (colorB.b - colorA.b);
    }

    else if (args->viz_type == 2){
      // Visaulizaing both //

      glm::vec4 colorA;
      glm::vec4 colorB;
    
      // Frequency visualization
      double freq = part->getFreq();
      
      assert( 20 <= freq && freq <= 20000); // Error if out of range

      double val = ( freq - 20 ) / ( 20000 - 20 ); // part/whole

      
      if( freq >= 10000){
      
        colorA = YELLOW;
        colorB = RED;
        val =  (freq-10000.0) / 10000.0;
      
      
      } else if(freq >= 1000){
      
        colorA = GREEN;
        colorB = YELLOW;
        val =  (freq-1000.0) / (10000-1000);
      
      } else if(freq >= 100){

        colorA = LITE_BLU;
        colorB = GREEN;
        val =  (freq-100.0) / (1000-100);

      } else {
      
        colorA = BLUE;
        colorB = LITE_BLU;
        val =  (freq-10.0) / (100-10);
      
      }
      
      color.r = colorA.r + val * (colorB.r - colorA.r);
      color.g = colorA.g + val * (colorB.g - colorA.g);
      color.b = colorA.b + val * (colorB.b - colorA.b);

      // Convert into dBs
      double dbs = 10 * log10( part->getWatt() / MIN_WATTAGE);
      double dbs_val = dbs / 120.0;

      // Sanity check
      if(dbs < 0 || dbs > 120){
        //std::cout << dbs << std::endl;
        dbs_val = 0;
      }
      
      color.a = dbs_val; // opacit
    }

    else if (args->viz_type == 3){

     // Visualizing wave fronts
     glm::vec3 dir = part->getDir();

     double rel_x =  (dir.x + 1) / 2.0;
     double rel_y =  (dir.y + 1) / 2.0;
     double rel_z =  (dir.z + 1) / 2.0;

     color = glm::vec4(rel_x, rel_y, rel_z,1);
    }

    else{

      // We will color by how best we fit to a shape
      // Colors we will use are
      // Not enough -- White
      // 4sides  = Red
      // 5sides  = Green
      // 6sides  = Blue
      // 7sides  = Magenta
      // 8+sides = Black

      // What I require to get value
      uint shape = getLowestCostShape(part);
      printf("Particle %p, Shape: %d\n",part,shape);
      switch(shape){
        case 0:
          color = glm::vec4(1,1,1,1); break; // white
        case 4:
          color = glm::vec4(1,0,0,1); break; // red
        case 5:
          color = glm::vec4(0,1,0,1); break; // green
        case 6:
          color = glm::vec4(0,0,1,1); break; // blue
        case 7:
          color = glm::vec4(0,0,0,1); break; // black
        default:
          color = glm::vec4(0,0,0,1); break; // black
      }

    }//vizType


    // DEBUG///////////////////////////////////////////////////////////////////
    // This code causes the animation to stop if we encounter a "free-particle"
    // a particle without having to hit any surface in the open.
    if(part->getCollisionSteps() == (int) pow(10,6) ){
      color = glm::vec4(1,0,0,1);
      std::cout << "Free Particle Caught in Code"<< std::endl;
      args->animate = false;
    }
    // END /////////////////////////////////////////////////////////////////////
    
    particle_verts.push_back(VBOPosNormalColor(pos,normal,color));
  
  }

  glBindBuffer(GL_ARRAY_BUFFER,particle_verts_VBO); 

  glBufferData(
      GL_ARRAY_BUFFER,
      sizeof(VBOPosNormalColor)*particle_verts.size(),
      &particle_verts[0],
      GL_DYNAMIC_DRAW); 
}

void ParticleSystem::drawParticles(){

  // Experimental Transparancy for Particles
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  HandleGLError("enter drawParticles");
  glPointSize(4) ;  // CHANGE ME <------------------------------- back to 2
  glBindBuffer(GL_ARRAY_BUFFER, particle_verts_VBO);

  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor),(void*)0);

  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );

  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 4, GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));

  glDrawArrays(GL_POINTS, 0, particle_verts.size());

  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);


  glDisable(GL_BLEND);

  HandleGLError("leaving drawParticles");
}


void ParticleSystem::setupCursorPoint(){

  glm::vec3 cursorPos = cursor;
  glm::vec3 normal(0,1,0);
  glm::vec4 color(1,0,0,1);

  // Binding the data
  cursor_verts.push_back( VBOPosNormalColor(cursorPos,normal,color));

  glBindBuffer(GL_ARRAY_BUFFER,cursor_verts_VBO); 

  glBufferData(
      GL_ARRAY_BUFFER,
      sizeof(VBOPosNormalColor)*1,
      &cursor_verts[0],
      GL_STATIC_DRAW); 
}

void ParticleSystem::drawCursorPoint(){

  HandleGLError("enter drawCursorPoint");
  glPointSize( 5 ); // <------------------- Added
  glBindBuffer(GL_ARRAY_BUFFER, cursor_verts_VBO);

  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor),(void*)0);

  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );

  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));

  glDrawArrays(GL_POINTS, 0, cursor_verts.size());

  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);

  HandleGLError("leaving drawCursorPoint");
}


void ParticleSystem::setupVelocityVisual(){

  HandleGLError("enter setupVelocityVisual");
  // Set up the points here
  float thickness = 0.001 * GLCanvas::bbox.maxDim();

  // For each particle
  for(int i = 0; i < particles.size(); i++){

    Particle * cur = particles[i];
    if(cur == NULL){
    }
    glm::vec3 pos = cur->getPos();
    glm::vec3 dir = cur->getDir();

    // std::cout << "================================\n";
    // std::cout << "Cur: " << *cur << std::endl;

    if(args->direction){
    
      addEdgeGeometry(velocity_verts,velocity_tri_indices,
          pos , // Start position
          pos +  ((float) 0.2) * dir, // Move a little in the dir
          glm::vec4(dir.x,dir.y,dir.z,1), 
          glm::vec4(dir.x,dir.y,dir.z,1), 
          thickness,
          thickness*0.1);
    }
  }

  // Create and setup velocity_tri_indices, velocity_verts
  // Setup the VBOS here
  glBindBuffer(GL_ARRAY_BUFFER,velocity_verts_VBO);
  glBufferData( 
      GL_ARRAY_BUFFER,
      sizeof(VBOPosNormalColor)*velocity_verts.size(),
      &velocity_verts[0],
      GL_STATIC_DRAW
      );

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,velocity_tri_indices_VBO);
  glBufferData(
      GL_ELEMENT_ARRAY_BUFFER,
      sizeof(VBOIndexedTri)*velocity_tri_indices.size(),
      &velocity_tri_indices[0],
      GL_STATIC_DRAW
      );

  HandleGLError("leave setupVelocityVisual");
}

void ParticleSystem::drawVelocityVisual(){

  HandleGLError("enter drawVelocityVisual");

  if(args->direction){
  
    glBindBuffer(GL_ARRAY_BUFFER, velocity_verts_VBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, velocity_tri_indices_VBO);
  
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE, 2*sizeof(glm::vec3) + sizeof(glm::vec4), (void*)0);

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE, 2*sizeof(glm::vec3) + sizeof(glm::vec4), (void*)sizeof(glm::vec3));

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2,3,GL_FLOAT,GL_FALSE, 2*sizeof(glm::vec3) + sizeof(glm::vec4), (void*)(sizeof(glm::vec3)*2));

    glDrawElements(GL_TRIANGLES, velocity_tri_indices.size()*3,
        GL_UNSIGNED_INT, BUFFER_OFFSET(0));

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);

  }

  HandleGLError("leave drawVelocityVisual");
}

void ParticleSystem::setupEdges(){
  // This function will setup the outline_vert vector of points

  assert(false); // No longer supported

  // If I have nothing to setup don't
  if(particles.size() == 0)
    return;

  int index = args->render_mask;

  // DEEBUG
  // if(index != 0 ) continue; // TODO REMOVE

  Particle * cur = particles[index];

  // GATHER STEP //////////////////////////////////////////////////////////

  float gather_distance = GATHER_DISTANCE;
  float gather_angle    = GATHER_ANGLE;

  std::vector<unsigned int> gathered_particles_indices;

  for(int i = 0; i < particles.size(); i++){

    Particle * other = particles[i];

    if( other == cur ) // dont count self
      continue;
    

    float dist = glm::distance(cur->getOldPos(), other->getOldPos());

    // Are we close enough
    if( dist < gather_distance){
    
      float angle = acos( glm::dot( cur->getDir(), other->getDir() ) / 
          (glm::length(cur->getDir()) * glm::length(other->getDir())));
      
      // Are we traveling together
      if( angle < gather_angle )
        gathered_particles_indices.push_back(i);
    
    }

  }

  // Find Mask Step ///////////////////////////////////////////////////////

  // Get particles that are alive for mask concideration
  std::vector < Particle * > particle_for_mask_calc;
  particle_for_mask_calc.push_back(cur);

  for( int i = 0; i < gathered_particles_indices.size(); i++)
   particle_for_mask_calc.push_back(particles[gathered_particles_indices[i]]);

  printf("Gathered %d Particles for Mask\n", particle_for_mask_calc.size());

  Mask mask;
  generateMask(particle_for_mask_calc, mask);

  mask.renderCost(happyness_verts);


  // As well as the other VBO
  HandleGLError("entering setupOutlineVisual");

  // Bind the happyness
  glBindBuffer(GL_ARRAY_BUFFER,happyness_verts_VBO); 

  glBufferData(
      GL_ARRAY_BUFFER,
      sizeof(VBOPosNormalColor)*happyness_verts.size(),
      &happyness_verts[0],
      GL_STATIC_DRAW); 

  HandleGLError("leaving setupOutlineVisual");
}

void ParticleSystem::drawHappinessVisual(){

  assert (false); // No longer supported
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  HandleGLError("enter drawHappyVisual");
  glBindBuffer(GL_ARRAY_BUFFER, happyness_verts_VBO);


  // glLineWidth(3);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor),(void*)0);

  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );

  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 4, GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));

  HandleGLError("enter debugDrawCall");
  glDrawArrays(GL_LINES, 0, happyness_verts.size());
  HandleGLError("leaving debugDrawCall");

  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);

  glDisable(GL_BLEND);

  // glLineWidth(3);
  HandleGLError("leaving drawHappyVisual");
}


void ParticleSystem::setupDelusionalParticles(){
  HandleGLError("entering setupDelusionalParticles");


  // If I have nothing to setup don't
  if(particles.size() == 0)
    return;

  // ==========================================================================
  // Gathering of the particles
  // ==========================================================================

  // Properties we will use for gathering 
  float gather_distance = GATHER_DISTANCE;
  float gather_angle    = GATHER_ANGLE;

  // Get the particle I will render for
  int index = args->render_mask;
  Particle * cur = particles[index];

  // We will will put particles we've gathered
  PartPtrVec gathered_particles;      
  gathered_particles.push_back(cur);

  // Use gather our particles
  particle_kdtree.GatherParticles(cur, 
    GATHER_DISTANCE, GATHER_ANGLE, gathered_particles);
  
  printf("gathered_particles.size() : %d \n", gathered_particles.size());

  // ==========================================================================
  // Setup of delusional particle points
  // ==========================================================================
  
  Mask mask;
  generateMask(gathered_particles, mask);

  int shape = mask.getShape();

  // Pushing in  center of mask
  delusional_verts.push_back(
      VBOPosNormalColor(cur->getOldPos(),glm::vec3(0,0,0),glm::vec4(1,1,1,1)));

  // Get those points
  std::vector < glm::vec3 > positions;
  printf("maskShape: %d\n", shape);
  delusionalParticleLocations(cur, gathered_particles,positions,shape);

  // ==========================================================================
  // Setup of the 1 -1 matching between delisional pos and our mask
  // ==========================================================================


  // get particles there
  const std::vector<Particle*> maskParticles = mask.getMaskParticles();

  // Assert 1-1 matching
  printf("maskParticles.size() = %d,  postions.size() = %d\n",maskParticles.size(), positions.size());
  assert(maskParticles.size()==positions.size());

  printf("Found Mask of Size: %d\n", mask.size() );

  for(int i = 0; i < maskParticles.size(); i++){

   Particle * curMask = maskParticles[i];

   if(curMask == NULL){
    delusional_verts.push_back(
        VBOPosNormalColor(positions[i],glm::vec3(0,0,0),glm::vec4(1,1,1,0.2)));
     continue;
   }
   
  glm::vec4 color = GiveHeapMapping( mask.getCost(i)/ (RADIUS_PARTICLE_WAVE * 1000.0));

  color.a = 0.1;
   connection_verts.push_back(
       VBOPosNormalColor(curMask->getOldPos(), glm::vec3(0,0,0), color));

  color.a = 1;
   connection_verts.push_back(
       VBOPosNormalColor(positions[i], glm::vec3(0,0,0), color));

  color.a = 0.5;
  delusional_verts.push_back(
      VBOPosNormalColor(positions[i],glm::vec3(0,0,0),color));
 }

  // Bind the delusional particles
  glBindBuffer(GL_ARRAY_BUFFER,delusional_verts_VBO); 

  
  glBufferData(
      GL_ARRAY_BUFFER,
      sizeof(VBOPosNormalColor)*delusional_verts.size(),
      &delusional_verts[0],
      GL_STATIC_DRAW); 

  // Bind the happyness
  glBindBuffer(GL_ARRAY_BUFFER,connection_verts_VBO); 

  glBufferData(
      GL_ARRAY_BUFFER,
      sizeof(VBOPosNormalColor)*connection_verts.size(),
      &connection_verts[0],
      GL_STATIC_DRAW); 

  // VBO SETUP
  HandleGLError("leaving setupDelusionalParticles");
}


void ParticleSystem::drawDelusionalParticles(){
  HandleGLError("entering drawDelusionalParticles");

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glPointSize(10) ; 
  glBindBuffer(GL_ARRAY_BUFFER, delusional_verts_VBO);


  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor),(void*)0);

  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );

  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 4, GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));


  glDrawArrays(GL_POINTS, 0, delusional_verts.size());

  HandleGLError("leaving debug");

  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);


  glDisable(GL_BLEND);


  HandleGLError("leaving drawDelusionalParticles");
}


void ParticleSystem::drawDelusionalConnections(){

  HandleGLError("entering drawDelusionalConnections");

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glBindBuffer(GL_ARRAY_BUFFER, connection_verts_VBO);

  // glLineWidth(3);
  glEnableVertexAttribArray(0);

  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor),(void*)0);

  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );

  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 4, GL_FLOAT,GL_FALSE,
      sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));

  glDrawArrays(GL_LINES, 0, connection_verts.size());

  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);

  glDisable(GL_BLEND);

  HandleGLError("leaving drawDelusionalConnections");
}

void ParticleSystem::setupSphere(){


  sphere.setup(10,10, sphere_verts,sphere_tri_indices );

  // Create and setup velocity_tri_indices, velocity_verts
  // Setup the VBOS here
  glBindBuffer(GL_ARRAY_BUFFER,sphere_verts_VBO);
  glBufferData( 
      GL_ARRAY_BUFFER,
      sizeof(VBOPosNormalColor)*sphere_verts.size(),
      &sphere_verts[0],
      GL_STATIC_DRAW );

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,sphere_tri_indices_VBO);

  glBufferData(
      GL_ELEMENT_ARRAY_BUFFER,
      sizeof(VBOIndexedTri)*sphere_tri_indices.size(),
      &sphere_tri_indices[0],
      GL_STATIC_DRAW );

  HandleGLError("leave setupVelocityVisual");

}

void ParticleSystem::drawSphere(){
  
  // This will setup the memeber functions so that you can render them
  // std::vector<VBOPosNormalColor> sphere_verts;
  // std::vector<VBOIndexedTri>     sphere_tri_indices;

  HandleGLError("enter drawSphereVisual");

  // if we are viewing our fitting types
  if(args->viz_type == 4){
  
    glBindBuffer(GL_ARRAY_BUFFER, sphere_verts_VBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphere_tri_indices_VBO);
  
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE, 2*sizeof(glm::vec3) + sizeof(glm::vec4), (void*)0);

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE, 2*sizeof(glm::vec3) + sizeof(glm::vec4), (void*)sizeof(glm::vec3));

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2,3,GL_FLOAT,GL_FALSE, 2*sizeof(glm::vec3) + sizeof(glm::vec4), (void*)(sizeof(glm::vec3)*2));

    glDrawElements(GL_TRIANGLES, sphere_tri_indices.size()*3,
        GL_UNSIGNED_INT, BUFFER_OFFSET(0));

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);

  }

  HandleGLError("leave drawSphereVisual");

}


















