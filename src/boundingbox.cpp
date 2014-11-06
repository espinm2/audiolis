#include <iostream>
#include "glCanvas.h"
#include "boundingbox.h"
#include "vbo_structs.h"
#include "render_utils.h"

// ====================================================================
// ====================================================================

void BoundingBox::initializeVBOs() {
  glGenBuffers(1, &bb_verts_VBO);
  glGenBuffers(1, &bb_tri_indices_VBO);
}


void BoundingBox::setupVBOs() {
  HandleGLError("bounding box setup VBOs enter");
  float thickness = 0.001*glm::length(maximum-minimum);

  glm::vec3& A = minimum;
  glm::vec3& B = maximum;
  glm::vec4 black(0,0,0,1);

  addEdgeGeometry(bb_verts,bb_tri_indices,glm::vec3(A.x,A.y,A.z),glm::vec3(A.x,A.y,B.z),black,black,thickness,thickness);
  addEdgeGeometry(bb_verts,bb_tri_indices,glm::vec3(A.x,A.y,B.z),glm::vec3(A.x,B.y,B.z),black,black,thickness,thickness);
  addEdgeGeometry(bb_verts,bb_tri_indices,glm::vec3(A.x,B.y,B.z),glm::vec3(A.x,B.y,A.z),black,black,thickness,thickness);
  addEdgeGeometry(bb_verts,bb_tri_indices,glm::vec3(A.x,B.y,A.z),glm::vec3(A.x,A.y,A.z),black,black,thickness,thickness);

  addEdgeGeometry(bb_verts,bb_tri_indices,glm::vec3(B.x,A.y,A.z),glm::vec3(B.x,A.y,B.z),black,black,thickness,thickness);
  addEdgeGeometry(bb_verts,bb_tri_indices,glm::vec3(B.x,A.y,B.z),glm::vec3(B.x,B.y,B.z),black,black,thickness,thickness);
  addEdgeGeometry(bb_verts,bb_tri_indices,glm::vec3(B.x,B.y,B.z),glm::vec3(B.x,B.y,A.z),black,black,thickness,thickness);
  addEdgeGeometry(bb_verts,bb_tri_indices,glm::vec3(B.x,B.y,A.z),glm::vec3(B.x,A.y,A.z),black,black,thickness,thickness);
  
  addEdgeGeometry(bb_verts,bb_tri_indices,glm::vec3(A.x,A.y,A.z),glm::vec3(B.x,A.y,A.z),black,black,thickness,thickness);
  addEdgeGeometry(bb_verts,bb_tri_indices,glm::vec3(A.x,A.y,B.z),glm::vec3(B.x,A.y,B.z),black,black,thickness,thickness);
  addEdgeGeometry(bb_verts,bb_tri_indices,glm::vec3(A.x,B.y,B.z),glm::vec3(B.x,B.y,B.z),black,black,thickness,thickness);
  addEdgeGeometry(bb_verts,bb_tri_indices,glm::vec3(A.x,B.y,A.z),glm::vec3(B.x,B.y,A.z),black,black,thickness,thickness);

  glBindBuffer(GL_ARRAY_BUFFER,bb_verts_VBO); 
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*bb_verts.size(),&bb_verts[0],GL_STATIC_DRAW); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,bb_tri_indices_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(VBOIndexedTri)*bb_tri_indices.size(),&bb_tri_indices[0],GL_STATIC_DRAW);

  HandleGLError("bounding box setup VBOs finished");
}

void BoundingBox::drawVBOs() {
  HandleGLError("draw VBOs a ");
  glBindBuffer(GL_ARRAY_BUFFER, bb_verts_VBO);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bb_tri_indices_VBO);
  glDrawElements(GL_TRIANGLES,
                 bb_tri_indices.size()*3,
                 GL_UNSIGNED_INT, 0);
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);
}

void BoundingBox::cleanupVBOs() {
  glDeleteBuffers(1, &bb_verts_VBO);
  glDeleteBuffers(1, &bb_tri_indices_VBO);
}


// ====================================================================
// ====================================================================
