#include "render_utils.h"
#include "glCanvas.h"


// since glLineWidth is gone...  
// instead we'll draw a rectangular box 
// (should probably use a geometry shader instead)
void addEdgeGeometry(std::vector<VBOPosNormalColor> &verts,
                     std::vector<VBOIndexedTri> &tri_indices,
                     const glm::vec3 &a, const glm::vec3 &b, 
                     const glm::vec4 &acolor, const glm::vec4 &bcolor, 
                     float a_th,float b_th) {
  
  // find perpendicular axes
  float length = glm::length(a-b);
  
  // std::cout << "length " << length << " vs " << std::min(a_th,b_th) << std::endl;
  if (length < 0.01*std::max(0.0000001f,std::min(a_th,b_th))) return;
  glm::vec3 dir = glm::normalize(b-a);
  glm::vec3 tmp = glm::cross(dir,glm::vec3(1,0,0));
  if (glm::length(tmp) < 0.1) {
    tmp = glm::cross(dir,glm::vec3(0,0,1));
  }
  tmp = glm::normalize(tmp);
  glm::vec3 one = glm::cross(dir,tmp);
  assert (fabs(glm::length(one)-1.0) < 0.001);
  glm::vec3 two = glm::cross(dir,one);
  assert (fabs(glm::length(two)-1.0) < 0.001);

  // draw the 6 faces of the box
  int start;
  start = verts.size();
  verts.push_back(VBOPosNormalColor(a-one*a_th+two*a_th,two,acolor));
  verts.push_back(VBOPosNormalColor(b-one*b_th+two*b_th,two,bcolor));
  verts.push_back(VBOPosNormalColor(b+one*b_th+two*b_th,two,bcolor));
  verts.push_back(VBOPosNormalColor(a+one*a_th+two*a_th,two,acolor));
  tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));
  tri_indices.push_back(VBOIndexedTri(start,start+2,start+3));
  start += 4;
  verts.push_back(VBOPosNormalColor(a-one*a_th-two*a_th,-two,acolor));
  verts.push_back(VBOPosNormalColor(b-one*b_th-two*b_th,-two,bcolor));
  verts.push_back(VBOPosNormalColor(b+one*b_th-two*b_th,-two,bcolor));
  verts.push_back(VBOPosNormalColor(a+one*a_th-two*a_th,-two,acolor));
  tri_indices.push_back(VBOIndexedTri(start,start+2,start+1));
  tri_indices.push_back(VBOIndexedTri(start,start+3,start+2));
  start += 4;
  verts.push_back(VBOPosNormalColor(a-two*a_th+one*a_th,one,acolor));
  verts.push_back(VBOPosNormalColor(b-two*b_th+one*b_th,one,bcolor));
  verts.push_back(VBOPosNormalColor(b+two*b_th+one*b_th,one,bcolor));
  verts.push_back(VBOPosNormalColor(a+two*a_th+one*a_th,one,acolor));
  tri_indices.push_back(VBOIndexedTri(start,start+2,start+1));
  tri_indices.push_back(VBOIndexedTri(start,start+3,start+2));
  start += 4;
  verts.push_back(VBOPosNormalColor(a-two*a_th-one*a_th,-one,acolor));
  verts.push_back(VBOPosNormalColor(b-two*b_th-one*b_th,-one,bcolor));
  verts.push_back(VBOPosNormalColor(b+two*b_th-one*b_th,-one,bcolor));
  verts.push_back(VBOPosNormalColor(a+two*a_th-one*a_th,-one,acolor));
  tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));
  tri_indices.push_back(VBOIndexedTri(start,start+2,start+3));
  start += 4;
  verts.push_back(VBOPosNormalColor(a-two*a_th-one*a_th,-dir,acolor));
  verts.push_back(VBOPosNormalColor(a-two*a_th+one*a_th,-dir,acolor));
  verts.push_back(VBOPosNormalColor(a+two*a_th+one*a_th,-dir,acolor));
  verts.push_back(VBOPosNormalColor(a+two*a_th-one*a_th,-dir,acolor));
  tri_indices.push_back(VBOIndexedTri(start,start+2,start+1));
  tri_indices.push_back(VBOIndexedTri(start,start+3,start+2));
  start += 4;
  verts.push_back(VBOPosNormalColor(b-two*b_th-one*b_th,dir,bcolor));
  verts.push_back(VBOPosNormalColor(b-two*b_th+one*b_th,dir,bcolor));
  verts.push_back(VBOPosNormalColor(b+two*b_th+one*b_th,dir,bcolor));
  verts.push_back(VBOPosNormalColor(b+two*b_th-one*b_th,dir,bcolor));
  tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));
  tri_indices.push_back(VBOIndexedTri(start,start+2,start+3));
}


#define HEAT_RED getColor(255,51,51,1)
#define HEAT_YELLOW getColor(255,255,51,1)
#define HEAT_GREEN getColor(51,255,51,1)
#define HEAT_BLUE getColor(51,51,255,1)
#define HEAT_BLACK getColor(0,0,0,1)

glm::vec4 GiveHeapMapping(double position){
  // assume it is [0,1]

  glm::vec4 colorA;
  glm::vec4 colorB;



  if(position > 0.80){
    // yellow -> red
    colorA = HEAT_YELLOW;
    colorB = HEAT_RED;
      
    

  }else if( position > 0.50){
    // green - > yellow
    colorA = HEAT_GREEN;
    colorB = HEAT_YELLOW;
  
  }else if (position > 0.1){
    // blue to green 
    colorA = HEAT_BLUE;
    colorB = HEAT_GREEN;
  
  
  }else{
    // black to blue
    colorA = HEAT_BLACK;
    colorB = HEAT_BLUE;
  
  }

  glm::vec4 color;
  color.r = colorA.r + (float)position * (colorB.r - colorA.r);
  color.g = colorA.g + (float)position * (colorB.g - colorA.g);
  color.b = colorA.b + (float)position * (colorB.b - colorA.b);

  std::cout << color.r <<  ", "  << color.g <<  ", "  << color.b  << std::endl;

  return color;

}



/* based on Delphi function by Witold J.Janik */
glm::vec4 GiveRainbowColor(double position)
{


  unsigned char c[] = {0,0,0};

  /* if position > 1 then we have repetition of colors it maybe useful    */
  if (position>1.0){if (position-(int)position==0.0)position=1.0; else position=position-(int)position;}
 
 
  unsigned char nmax=6; /* number of color segments */
  double m=nmax* position;
 
  int n=(int)m; // integer of m
 
  double f=m-n;  // fraction of m
  unsigned char t=(int)(f*255);
 
  switch( n){
     case 0: {
        c[0] = 255;
        c[1] = t;
        c[2] = 0;
         break;
      };
     case 1: {
        c[0] = 255 - t;
        c[1] = 255;
        c[2] = 0;
         break;
      };
     case 2: {
        c[0] = 0;
        c[1] = 255;
        c[2] = t;
         break;
      };
     case 3: {
        c[0] = 0;
        c[1] = 255 - t;
        c[2] = 255;
         break;
      };
     case 4: {
        c[0] = t;
        c[1] = 0;
        c[2] = 255;
         break;
      };
     case 5: {
        c[0] = 255;
        c[1] = 0;
        c[2] = 255 - t;
         break;
      };
      default: {
        c[0] = 255;
        c[1] = 0;
        c[2] = 0;
         break;
      };
  }; // case



  return glm::vec4(c[0]/256.0 , c[1] / 256.0, c[2] /256.0, 1);
}

