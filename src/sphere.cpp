#include "sphere.h"
#include <glm/glm.hpp>

// helper function to place a grid of points on the sphere
glm::vec3 Sphere::ComputeSpherePoint(float s, float t, const glm::vec3 center, float radius) {
  float angle = 2*M_PI*s;
  float y = -cos(M_PI*t);
  float factor = sqrt(1-y*y);
  float x = factor*cos(angle);
  float z = factor*-sin(angle);
  glm::vec3 answer = glm::vec3(x,y,z);
  answer *= radius;
  answer += center;
  return answer;
}

void Sphere::setup(
  int sphere_horiz, int sphere_vert,
  std::vector<VBOPosNormalColor> & verts, 
  std::vector<VBOIndexedTri> & tri_indices ){
  // inputs
  //      verts | where our verticies are being stored
  //    tri_indices |  where we store our triangle indices
  //
  // outputs
  //      All is returned by refrence
  
  // and convert it into quad patches for radiosity
  int h = sphere_horiz; int v = sphere_vert;
  assert (h % 2 == 0);
  int i,j,start; int va,vb,vc,vd;
  glm::vec3 a,b,c,d; std::vector<glm::vec3> spherePoints;
  glm::vec4 color(0.3,0.3,0.3,1);

  // place vertices into my temp array spherePoints
  spherePoints.push_back(center+radius*glm::vec3(0,-1,0));  // bottom
  for (j = 1; j < v; j++) {  // middle
    for (i = 0; i < h; i++) {
      float s = i / float(h);
      float t = j / float(v);
      spherePoints.push_back(ComputeSpherePoint(s,t,center,radius));
    }
  }

  spherePoints.push_back(center+radius*glm::vec3(0,1,0));  // top

  // the middle patches
  for (j = 1; j < v-1; j++) {
    for (i = 0; i < h; i++) {

      va = 1 +  i      + h*(j-1);
      vb = 1 + (i+1)%h + h*(j-1);
      vc = 1 +  i      + h*(j);
      vd = 1 + (i+1)%h + h*(j);

      a = spherePoints[va];
      b = spherePoints[vb];
      c = spherePoints[vc];
      d = spherePoints[vd];


      // Create the two resulting triangles
      start = verts.size();
      verts.push_back(VBOPosNormalColor(a, glm::normalize(a-center) , color ) );
      verts.push_back(VBOPosNormalColor(b, glm::normalize(b-center) , color ) );
      verts.push_back(VBOPosNormalColor(c, glm::normalize(c-center) , color ) );
      verts.push_back(VBOPosNormalColor(d, glm::normalize(d-center) , color ) );
      tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));
      tri_indices.push_back(VBOIndexedTri(start+1,start+2,start+3));

    }
  }

  for (i = 0; i < h; i+=2) {
    
    // the bottom patches
    va = 0;
    vb = 1 +  i;
    vc = 1 + (i+1)%h;
    vd = 1 + (i+2)%h;

    a = spherePoints[va];
    b = spherePoints[vb];
    c = spherePoints[vc];
    d = spherePoints[vd];

    start = verts.size();
    verts.push_back(VBOPosNormalColor(d, glm::normalize(a-center) , color ) );
    verts.push_back(VBOPosNormalColor(c, glm::normalize(b-center) , color ) );
    verts.push_back(VBOPosNormalColor(b, glm::normalize(c-center) , color ) );
    verts.push_back(VBOPosNormalColor(a, glm::normalize(d-center) , color ) );
    tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));
    tri_indices.push_back(VBOIndexedTri(start,start+2,start+3));

    // the top patches
    va = 1 + h*(v-1);
    vb = 1 +  i      + h*(v-2);
    vc = 1 + (i+1)%h + h*(v-2);
    vd = 1 + (i+2)%h + h*(v-2);

    a = spherePoints[va];
    b = spherePoints[vb];
    c = spherePoints[vc];
    d = spherePoints[vd];

    start = verts.size();
    verts.push_back(VBOPosNormalColor(b, glm::normalize(a-center) , color ) );
    verts.push_back(VBOPosNormalColor(c, glm::normalize(b-center) , color ) );
    verts.push_back(VBOPosNormalColor(d, glm::normalize(c-center) , color ) );
    verts.push_back(VBOPosNormalColor(a, glm::normalize(d-center) , color ) );
    tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));
    tri_indices.push_back(VBOIndexedTri(start,start+2,start+3));

  }
}
