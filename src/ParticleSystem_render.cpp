#include <glm/glm.hpp>

#include "ParticleSystem.h"
#include "particle.h"
#include "boundingbox.h"
#include "camera.h"

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
    glm::vec3 normal = part->getDir(); 
    
    
    // Picking color
    glm::vec4 color(0,0,1,1);
    
    particle_verts.push_back(VBOPosNormalColor(pos,normal,color));
  
  }

  glBindBuffer(GL_ARRAY_BUFFER,particle_verts_VBO); 

  glBufferData(
      GL_ARRAY_BUFFER,
      sizeof(VBOPosNormalColor)*particle_verts.size(),
      &particle_verts[0],
      GL_STATIC_DRAW); 
}

void ParticleSystem::drawParticles(){


  HandleGLError("enter drawParticles");
  glPointSize( 2 ) ; 
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

  HandleGLError("leaving drawParticles");

}


void ParticleSystem::setupCursorPoint(){

  glm::vec3 cursorPos( cursor.x(), cursor.y(), cursor.z() );
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
  glPointSize( getGLPointSize(cursor) ); // <------------------- Added
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


int ParticleSystem::getGLPointSize(const glm::vec3 & point){

  // Get Camera  Position
  glm::vec3 cameraPos = GLCanvas::camera->camera_position;

  double dist = glm::distance(point,cameraPos);

  return 10;


}
