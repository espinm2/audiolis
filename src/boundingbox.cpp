#include <iostream>
#include "glCanvas.h"
#include "boundingbox.h"
#include "vbo_structs.h"
#include "render_utils.h"
#include "ray.h"

// ====================================================================
// ====================================================================

// Given a ray will return true or false if bbox was hit
bool BoundingBox::hitbox(const Ray &r, double t1, double t2){

  // Check to see if this ray hits me
  glm::vec3 d = r.getDirection();
  glm::vec3 e = r.getOrigin();
  double velocity = r.getVelocity();

  // Intersection calculations
  double tx_min, tx_max, // Assign of these
         ty_min, ty_max,
         tz_min, tz_max,
         a;

  // X's
  a = 1.0 / (d.x * velocity) ;
  if(a >= 0 ){

    tx_min =  a * ( getMin().x - e.x );
    tx_max =  a * ( getMax().x - e.x );

  }else{

    tx_min =  a * ( getMax().x - e.x );
    tx_max =  a * ( getMin().x - e.x );
  
  }

  // Y's
  a = 1.0 / ( d.y * velocity );
  if(a >= 0 ){

    ty_min =  a * ( getMin().y - e.y );
    ty_max =  a * ( getMax().y - e.y );

  }else{

    ty_min =  a * ( getMax().y - e.y );
    ty_max =  a * ( getMin().y - e.y );
  
  }


  // Z's
  a = 1.0 / ( d.z * velocity );
  if(a >= 0 ){

    tz_min =  a * ( getMin().z - e.z );
    tz_max =  a * ( getMax().z - e.z );

  }else{

    tz_min =  a * ( getMax().z - e.z );
    tz_max =  a * ( getMin().z - e.z );
  
  }

  // Check for overlap
  bool x_y_inter, y_z_inter, z_x_inter;
  x_y_inter = tx_min <= ty_max  && ty_min <= tx_max;
  y_z_inter = ty_min <= tz_max  && tz_min <= ty_max;
  z_x_inter = tz_min <= tx_max  && tx_min <= tz_max;

  // bool x_sat = tx_min <= time_step && time_step <= tx_max;
  // bool y_sat = ty_min <= time_step && time_step <= ty_max;
  // bool z_sat = tz_min <= time_step && time_step <= tz_max;

  // If we do not share a common interval, no collision ever
  if(!x_y_inter || !y_z_inter || !z_x_inter){
    return false;
  }

  // find exactly when we hit our box
  double t_in = std::max( std::max( tx_min, ty_min), tz_min );
  double t_out = std::min( std::min( tx_max, ty_max), tz_max );



  // Assert this ordering
  assert(t_in <= t_out);

  // Do we hit in our timestep?
  bool collision = t_in <= t2  &&  t1 <= t_out;

  return collision;
    
}//end



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
