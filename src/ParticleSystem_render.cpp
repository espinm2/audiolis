#include <glm/glm.hpp>
#include <math.h>

#include "ParticleSystem.h"
#include "particle.h"
#include "boundingbox.h"
#include "camera.h"
#include "argparser.h"
#include "render_utils.h"

#define MAX_ITERATIONS 6000
#define MAX_ITERATIONS 6000

#define MIN_WATTAGE 0.000000000002

/*
// This is for rainbow color
#define BLUE     getColor(1,66,255,1)
#define LITE_BLU getColor(0,254,252,1)
#define GREEN    getColor(135,255,119,1)
#define YELLOW   getColor(255,240,1,1)
#define RED      getColor(255,30,0,1)
*/

#define BLUE     getColor(255,255,178,1)
#define LITE_BLU getColor(254,204,92,1)
#define GREEN    getColor(253,141,60,1)
#define YELLOW   getColor(240,59,31,1)
#define RED      getColor(189,0,38,1)


void ParticleSystem::initializeVBOs(){

  // Get uniquie ids from buffers
  glGenBuffers(1,&particle_verts_VBO);
  glGenBuffers(1,&cursor_verts_VBO);

}

void ParticleSystem::setupVBOs(){

  // Delete old data
  cursor_verts.clear();
  particle_verts.clear();


  // Setup new Data
  setupCursorPoint();
  setupParticles(); 

}

void ParticleSystem::drawVBOs(){
  HandleGLError("enter draw vbos");


  drawCursorPoint();
  drawParticles();

  HandleGLError("leaving draw vbos");

}

void ParticleSystem::cleanupVBOs(){
  glDeleteBuffers(1,&particle_verts_VBO);
  glDeleteBuffers(1,&cursor_verts_VBO);

}

void ParticleSystem::setupParticles(){

  // Push everything into particle_verts
  for(int i = 0 ; i < particles.size(); i++){

    // Cur Particle
    Particle * part = particles[i];

    // Getting pos
    glm::vec3 pos = part->getPos();

    // Picking Normal
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
      // Visaulizaing both ////////////////////////////////////////////////////

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
      
      color.a = dbs_val; // opacity
      // color.a = 1.0; // opacity // TEST

    }else{
    
      // Should error out
      std::cout << "Bad Error Type" << std::endl;
      assert(false);
    
    }

    if(part->getMaterialHit() == "none"){

      color = getColor(0,0,0,0);
    
    }
    
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
  glPointSize( 10 ); // <------------------- Added
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


