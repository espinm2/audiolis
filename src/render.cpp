#include <iostream>

#include "glCanvas.h"
#include "mesh.h"
#include "edge.h"
#include "vertex.h"
#include "triangle.h"
#include "argparser.h"
#include "utils.h"
#include <list>
#include <algorithm>
#include <math.h>       /* fabs */
#include <stdio.h>      /* printf */

// Predefined colors to use
glm::vec4 floor_color(0.9,0.8,0.7,1);
glm::vec4 mesh_color(0.8,0.8,0.8,1);
glm::vec4 mirror_color(0.1,0.1,0.2,1);
glm::vec4 mirror_tint(0.85,0.9,0.95,1);

glm::vec4 red(1.0,0,0,1);
glm::vec4 green(0,1,0,0.5);

float floor_factor = 0.75;

// =======================================================================
// =======================================================================


// the light position can be animated
glm::vec3 Mesh::LightPosition() const {
  glm::vec3 min = bbox.getMin();
  glm::vec3 max = bbox.getMax();
  glm::vec3 tmp;
  bbox.getCenter(tmp);
  tmp += glm::vec3(0,5.0*(max.y-min.y),0);
  tmp += glm::vec3(cos(args->timer) * (max.x-min.x), 0, 0);
  tmp += glm::vec3(0,0,sin(args->timer) * (max.z-min.z));
  // tmp = GLCanvas::camera->camera_position;
  return tmp;
}


void Mesh::initializeVBOs() {

  // Regular mesh buffer
  glGenBuffers(1,&mesh_tri_verts_VBO);
  glGenBuffers(1,&mesh_tri_indices_VBO);


  // Light buffer
  glGenBuffers(1,&light_vert_VBO);
  bbox.initializeVBOs();
}

void Mesh::cleanupVBOs() {
  glDeleteBuffers(1,&mesh_tri_verts_VBO);
  glDeleteBuffers(1,&mesh_tri_indices_VBO);
  glDeleteBuffers(1,&light_vert_VBO);
  bbox.cleanupVBOs();
}

// ================================================================================
// ================================================================================

void Mesh::SetupLight(const glm::vec3 &light_position) {
  light_vert.push_back(VBOPosNormalColor(light_position,glm::vec3(1,0,0),glm::vec4(1,1,0,0)));
  glBindBuffer(GL_ARRAY_BUFFER,light_vert_VBO); 
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*1,&light_vert[0],GL_STATIC_DRAW); 
}


void Mesh::SetupMesh() {
  for (triangleshashtype::iterator iter = triangles.begin();
       iter != triangles.end(); iter++) {
    Triangle *t = iter->second;
    glm::vec3 a = (*t)[0]->getPos();
    glm::vec3 b = (*t)[1]->getPos();
    glm::vec3 c = (*t)[2]->getPos();    
    glm::vec3 na = ComputeNormal(a,b,c);
    glm::vec3 nb = na;
    glm::vec3 nc = na;


    if( -1 <= na.y && na.y <= -1 + 0.01){
      continue;
    }


    if (args->gouraud_normals) {
      na = (*t)[0]->getGouraudNormal();
      nb = (*t)[1]->getGouraudNormal();
      nc = (*t)[2]->getGouraudNormal();
    }

    glm::vec4 color;
    color = glm::vec4(mesh_color.r * 0.1, mesh_color.g * 0.1, mesh_color.b*0.1,1);

    TriVBOHelper(mesh_tri_verts,mesh_tri_indices,
                 a,b,c,
                 na,nb,nc,
                 color,color,color
                 );

    /*
    int start = mesh_tri_verts.size();
    mesh_tri_verts.push_back(VBOPosNormalColor(a,na,mesh_color));
    mesh_tri_verts.push_back(VBOPosNormalColor(b,nb,mesh_color));
    mesh_tri_verts.push_back(VBOPosNormalColor(c,nc,mesh_color));
    mesh_tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));
    */
  }

  glBindBuffer(GL_ARRAY_BUFFER,mesh_tri_verts_VBO); 
  glBufferData(GL_ARRAY_BUFFER,
	       sizeof(VBOPosNormalColor) * mesh_tri_verts.size(), 
	       &mesh_tri_verts[0],
	       GL_STATIC_DRAW); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_tri_indices_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
	       sizeof(VBOIndexedTri) * mesh_tri_indices.size(),
	       &mesh_tri_indices[0], GL_STATIC_DRAW);


}



// ================================================================================
// ================================================================================

void Mesh::DrawLight() {
  HandleGLError("enter draw mirror");
  glPointSize(10);
  glBindBuffer(GL_ARRAY_BUFFER, light_vert_VBO);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
  glDrawArrays(GL_POINTS, 0, light_vert.size());
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);
  HandleGLError("enter draw mirror");
}


void Mesh::DrawMesh() {
  HandleGLError("enter draw mesh");
  glBindBuffer(GL_ARRAY_BUFFER,mesh_tri_verts_VBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_tri_indices_VBO); 

  // This is the posititons
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,
                        3,
                        GL_FLOAT,
                        GL_FALSE,
                        sizeof(VBOPosNormalColor),
                        (void*)0
                        );

  // This is the normal
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,
                        3,
                        GL_FLOAT,
                        GL_FALSE,
                        sizeof(VBOPosNormalColor),
                        (void*)sizeof(glm::vec3)
                        );

  // Is this the color?
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2,
                        3,
                        GL_FLOAT,
                        GL_FALSE,
                        sizeof(VBOPosNormalColor),
                        (void*)(sizeof(glm::vec3)*2)
                        );
  glDrawElements(GL_TRIANGLES, mesh_tri_indices.size()*3,GL_UNSIGNED_INT, 0);


  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);
  HandleGLError("leaving draw mesh");
}

// ======================================================================================
// ======================================================================================

void Mesh::setupVBOs() {
  // delete all the old geometry
  mesh_tri_verts.clear(); 
  mesh_tri_indices.clear();
  light_vert.clear();

  //edit
  extend_edges.clear();

  // setup the new geometry
  glm::vec3 light_position = LightPosition();
  SetupLight(light_position);
  SetupMesh();
  bbox.setupVBOs();

}

void Mesh::drawVBOs() {


  // mode 1: STANDARD PHONG LIGHTING (LIGHT ON)
  glUniform1i(GLCanvas::colormodeID, 1);

  // shader 0: NO SHADER
  glUniform1i(GLCanvas::whichshaderID, 0);


  HandleGLError("enter draw vbos");

  // --------------------------
  glDisable(GL_STENCIL_TEST);
  glClearDepth(1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);
  if (args->geometry) {
    // shader 1: CHECKERBOARD
    // shader 2: ORANGE
    // shader 3: other
    //
    glUniform1i(GLCanvas::whichshaderID, args->whichshader);
    DrawMesh();
    glUniform1i(GLCanvas::whichshaderID, 0);
  }


  // mode 0: NO LIGHTING
  glUniform1i(GLCanvas::colormodeID, 0);

  DrawLight();
  
  if (args->bounding_box) {
    bbox.drawVBOs();
  }

  HandleGLError(); 

}

// =================================================================
void Mesh::TriVBOHelper(std::vector<VBOPosNormalColor> &mesh_tri_verts,
                     std::vector<VBOIndexedTri> &mesh_tri_indices,
                     const glm::vec3 &pos_a,
                     const glm::vec3 &pos_b,
                     const glm::vec3 &pos_c,
                     const glm::vec3 &normal_a,
                     const glm::vec3 &normal_b,
                     const glm::vec3 &normal_c,
                     const glm::vec4 &color_ab,
                     const glm::vec4 &color_bc,
                     const glm::vec4 &color_ca)
{
    // To create a wireframe rendering...
  // Each mesh triangle is actually rendered as 3 small triangles
  //           b
  //          /|\
  //         / | \
  //        /  |  \
  //       /   |   \
  //      /    |    \
  //     /    .'.    \
  //    /  .'     '.  \
  //   /.'           '.\
  //  a-----------------c
  //

  // Final all the stuff for the centriod


  int start;

  if (args->wireframe) {

    glm::vec4 center_color(1,1,1,1);
    glm::vec3 centroid = 1.0f / 3.0f * (pos_a + pos_b + pos_c);
    glm::vec3 normal = normal_a + normal_b + normal_c;
    glm::normalize(normal);

    // WIREFRAME

    // make the 3 small triangles
    start = mesh_tri_verts.size();
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_a,normal_a,color_ab));
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_b,normal_b,color_ab));
    mesh_tri_verts.push_back(VBOPosNormalColor(centroid,normal,center_color));
    mesh_tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));


    // make the 3 small triangles
    start = mesh_tri_verts.size();
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_b,normal_b,color_bc));
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_c,normal_c,color_bc));
    mesh_tri_verts.push_back(VBOPosNormalColor(centroid,normal,center_color));
    mesh_tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));



    // make the 3 small triangles
    start = mesh_tri_verts.size();
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_c,normal_c,color_ca));
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_a,normal_a,color_ca));
    mesh_tri_verts.push_back(VBOPosNormalColor(centroid,normal,center_color));
    mesh_tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));


  } else {
    // NON WIREFRAME

    start = mesh_tri_verts.size();
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_a,normal_a,mesh_color));
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_b,normal_b,mesh_color));
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_c,normal_c,mesh_color));

    mesh_tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));

  }


}
