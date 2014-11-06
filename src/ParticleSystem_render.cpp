#include "ParticleSystem.h"
#include "boundingbox.h"
#include "camera.h"

void ParticleSystem::initializeVBOs(){

  // Get uniquie ids from buffers
  glGenBuffers(1,&particle_verts_VBO);
  glGenBuffers(1,&cursor_verts_VBO);

}

void ParticleSystem::setupVBOs(){

  // Delete old data
  particle_verts.clear();
  cursor_verts.clear();


  // Setup new Data
  // setupParticles(); TODO
  setupCursorPoint();
}

void ParticleSystem::drawVBOs(){
  //TODO Implement
  HandleGLError("enter draw vbos");

  /*
  glDisable(GL_STENCIL_TEST);
  glClearDepth(1);

  glClear(GL_COLOR_BUFFER_BIT 
  |GL_DEPTH_BUFFER_BIT | 

  GL_STENCIL_BUFFER_BIT);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);
  */

  // drawParticles(); TODO
  drawCursorPoint();

  HandleGLError("leaving draw vbos");

}

void ParticleSystem::cleanupVBOs(){
  //TODO Implement
  glDeleteBuffers(1,&particle_verts_VBO);
  glDeleteBuffers(1,&cursor_verts_VBO);

}

void ParticleSystem::setupParticles(){
  //TODO Implement
}

void ParticleSystem::drawParticles(){
  //TODO Implement
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


int ParticleSystem::getGLPointSize(const Vec3f & point){

  // Get Camera  Position
  glm::vec3 cPos = GLCanvas::camera->camera_position;
  Vec3f cameraPos(cPos.x, cPos.y, cPos.z);

  double dist = point.Distance3f(cameraPos);
  std::cout << dist << std::endl;


  return 10;


}
