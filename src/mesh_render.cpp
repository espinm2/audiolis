#include <iostream>

#include "glCanvas.h"
#include "mesh.h"
#include "edge.h"
#include "vertex.h"
#include "triangle.h"
#include "argparser.h"
#include "render_utils.h"
#include "camera.h"
#include "MersenneTwister.h"
#include "geometry_utils.h"
#include <list>
#include <algorithm>
#include <math.h>       /* fabs */
#include <stdio.h>      /* printf */
#include <cstdlib> 

// Predefined colors to use
glm::vec4 mesh_color(0.8,0.8,0.8,1);


// =======================================================================
// =======================================================================


// the light position can be animated
glm::vec3 Mesh::LightPosition() const {

  // glm::vec3 min = bbox.getMin();
  // glm::vec3 max = bbox.getMax();
  // glm::vec3 tmp;
  // bbox.getCenter(tmp);
  // tmp += glm::vec3(0,5.0*(max.y-min.y),0);
  // tmp += glm::vec3(cos(args->timer) * (max.x-min.x), 0, 0);
  // tmp += glm::vec3(0,0,sin(args->timer) * (max.z-min.z));

  // Orginal, hover camera over top of the scene
  glm::vec3 min = bbox.getMin();
  glm::vec3 max = bbox.getMax();
  float  radius = std::max((max.x-min.x), (max.z-min.z));

  // Get center
  glm::vec3 t;
  bbox.getCenter(t);

  // Get camera dir from center without y comp
  glm::vec3 dir = GLCanvas::camera->camera_position - t;
  dir = glm::normalize(glm::vec3(dir.x,0,dir.z));

  // Add offset
  // Fing projected
  t = t + (  dir * radius );
  
  // Add Y offset
  t += glm::vec3(0,5.0*(max.y-min.y),0);


  // glm::vec3 tmp = GLCanvas::camera->camera_position;
  return t;
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

  light_vert.push_back(
      VBOPosNormalColor(light_position,glm::vec3(1,0,0),glm::vec4(1,1,0,0)));

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
    glm::vec3 na = compute_normal(a,b,c);
    glm::vec3 nb = na;
    glm::vec3 nc = na;


    // TODO Make Render Floor Option
    
    if(args->render_top)
      if( -1 <= na.y && na.y <= -1 + 0.001)
        continue;

    // Looks really bad yo
    if (args->gouraud_normals) {
      na = (*t)[0]->getGouraudNormal();
      nb = (*t)[1]->getGouraudNormal();
      nc = (*t)[2]->getGouraudNormal();
    }

    // Color of triangle
    glm::vec4 center_color;
    glm::vec4 wire_color;

    std::string mtl = t->getMaterial();
    
    //unsigned int WALL_MATERIAL = 0; // This is a brick wall
    unsigned int WALL_MATERIAL = 1; // This is a Concrete wall
    // unsigned int WALL_MATERIAL = 2; // Ceramnic-Tiled wall

    unsigned int FLOOR_MATERIAL = 0; // pvc floor
     //unsigned int FLOOR_MATERIAL = 1; // carpeted floor

    if( mtl.compare(0,5,"GLASS") == 0){
      center_color = getColor(141,211,199,1);
    }

    // We will assume only some walls is made of Parete Absorber
    else if(mtl.compare(0,4,"wall") == 0){

        // Get the index
        std::string temp  = mtl.substr(5);
        int index = atoi(temp.c_str()) % 3;

        if(index == 1) {

          // Make this wall an absorber
         center_color = getColor(179,222,105,1);

        }else{

          // Make these walls anything
          switch (WALL_MATERIAL) {
            case 0:
              center_color = getColor(251,128,114,1);
              break;
            case 1:
              center_color = getColor(190,186,218,1);
              break;
            case 2:
              center_color = mesh_color;
              break;
            default:
              assert(false);
          }

        }
    } // walls

    else if( mtl.compare(0,6,"FILLIN") == 0) {

      switch (WALL_MATERIAL) {
        case 0:
          center_color = getColor(251,128,114,1);
          break;
        case 1:
          center_color = getColor(190,186,218,1);
          break;
        case 2:
          center_color = mesh_color;
          break;
        default:
          assert(false);
      }

    }

    else if( mtl == "floor"){
    
      switch (FLOOR_MATERIAL) {
        case 0:
          center_color = getColor(255,255,179,1);
          break;
        case 1:
          center_color = getColor(128,177,211,1);
          break;
        default:
          assert(false);
      }
    
    }

    else{
      // Assume ceiling
      // center_color = getColor(255,255,51,1);
      center_color = mesh_color;
    }
    // Finding colors OLD VERSION

    // if(mtl.compare(0,4,"wall") == 0){

    //   // Get the index
    //   std::string temp  = mtl.substr(5);

    //   int colorIndex = atoi(temp.c_str()) % 8;

    //   center_color = getColor(100,100,100,1);

    //   // DEBUG LATER (Colors wont change)
    //   
    //   switch (colorIndex) {
    //     case 0:
    //       center_color = getColor(228,26,28,1);
    //       break;
    //     case 1:
    //       center_color = getColor(247,129,191,1);
    //       break;
    //     case 2:
    //       center_color = getColor(55,126,184,1);
    //       break;
    //     case 3:
    //       center_color = getColor(77,175,74,1);
    //       break;
    //     case 4:
    //       center_color = getColor(152,78,163,1);
    //       break;
    //     case 5:
    //       center_color = getColor(255,127,0,1);
    //       break;
    //     case 6:
    //       center_color = getColor(255,255,51,1);
    //       break;
    //     case 7:
    //       center_color = getColor(166,86,40,1);
    //       break;
    //     default:
    //       assert(false);
    //   }

    //   wire_color = center_color;
    // 
    // }else if(mtl.compare(0,5,"GLASS") == 0){ 

    //   center_color = glm::vec4(0.00,1.00,0.80,1);
    //   wire_color = center_color;
    //   
    // }else if(mtl.compare(0,6,"FILLIN") == 0){ 

    //   center_color = glm::vec4(0.2,0.3,0,1);
    //   wire_color = center_color;
    //   
    // }else{
    //   
    //   center_color = mesh_color;

    //   wire_color = glm::vec4(
    //       mesh_color.r * 0.1, 
    //       mesh_color.g * 0.1, 
    //       mesh_color.b*0.1,
    //       1);

    // }

    wire_color = center_color; // <------------------------------------------------- Recall we killed the wiremesh

    // Sending color
    TriVBOHelper(mesh_tri_verts,mesh_tri_indices,
                 a,b,c,
                 na,nb,nc,
                 center_color,
                 wire_color,wire_color,wire_color);

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


  HandleGLError("enter draw vbos");

  // --------------------------
  glDisable(GL_STENCIL_TEST);
  glClearDepth(1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);
  if (args->geometry) {
    DrawMesh();
  }

  // If you want to draw light as glPoint
  // Keeping as an example to use glPoint
  // DrawLight();
  
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
                     const glm::vec4 &center_color,
                     const glm::vec4 &color_ab,
                     const glm::vec4 &color_bc,
                     const glm::vec4 &color_ca) {


  // Center color gets rendered in non-wireframe visualization
  // Otherwise we render the outlines provided
  
  int start;

  if (args->wireframe) {

    glm::vec3 centroid = 1.0f / 3.0f * (pos_a + pos_b + pos_c);
    glm::vec3 normal = normal_a + normal_b + normal_c;
    glm::normalize(normal);

    // WIREFRAME

    // make the 3 small triangles
    start = mesh_tri_verts.size();
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_a,normal_a,color_ab));
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_b,normal_b,color_ab));
    mesh_tri_verts.push_back(VBOPosNormalColor(centroid,normal,mesh_color));
    mesh_tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));


    // make the 3 small triangles
    start = mesh_tri_verts.size();
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_b,normal_b,color_bc));
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_c,normal_c,color_bc));
    mesh_tri_verts.push_back(VBOPosNormalColor(centroid,normal,mesh_color));
    mesh_tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));



    // make the 3 small triangles
    start = mesh_tri_verts.size();
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_c,normal_c,color_ca));
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_a,normal_a,color_ca));
    mesh_tri_verts.push_back(VBOPosNormalColor(centroid,normal,mesh_color));
    mesh_tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));


  } else {

    // NON WIREFRAME

    start = mesh_tri_verts.size();
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_a,normal_a,center_color));
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_b,normal_b,center_color));
    mesh_tri_verts.push_back(VBOPosNormalColor(pos_c,normal_c,center_color));

    mesh_tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));

  }


}
