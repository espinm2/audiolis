#ifndef _RENDER_UTILS_H
#define _RENDER_UTILS_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <vector>

#include "vbo_structs.h"


// =========================================================================
// These two functions convert between linear intensity values
// (approximate range 0->1) to an sRGB value (approximate range 0->1).
// The sRGB values make best use of 8 bit storage and are the best
// input for most displays and darkened viewing environments.

#define SRGB_ALPHA 0.055
#define BUFFER_OFFSET(i) ((char *)NULL + (i))

inline float linear_to_srgb(float x) {
  float answer;
  if (x <= 0.0031308)
    answer = 12.92*x;
  else 
    answer = (1+SRGB_ALPHA)*(pow(x,1/2.4)-SRGB_ALPHA);
  return answer;
}

inline float srgb_to_linear(float x) {
  float answer;
  if (x <= 0.04045)
    answer = x/12.92;
  else 
    answer = pow((x+SRGB_ALPHA)/(1+SRGB_ALPHA),2.4);
  return answer;
}

// =========================================================================
// utility functions 
inline glm::vec4 getColor(int r, int g, int b, float a){
  return glm::vec4(r/255.0, g/255.0, b/255.0, a);
}

// =========================================================================
// for rendering "lines"

void addEdgeGeometry(std::vector<VBOPosNormalColor> &verts,
                     std::vector<VBOIndexedTri> &tri_indices,
                     const glm::vec3 &a, const glm::vec3 &b, 
                     const glm::vec4 &acolor, const glm::vec4 &bcolor, 
                     float a_th,float b_th);

// =========================================================================
// helper functions for printing & reading vec3's

inline std::ostream& operator<<(std::ostream& ostr, const glm::vec3 &v) {
  ostr << "<" << v.x << "," << v.y << "," << v.z << ">";
  return ostr;
}

inline std::ostream& operator<<(std::ostream& ostr, const glm::vec4 &v) {
  ostr << "<" << v.r << "," << v.g << "," << v.b << ","  <<  v.a << ">";
  return ostr;
}

inline std::istream& operator>>(std::istream& istr, glm::vec3 &v) {
  char c;
  istr >> c;  assert (c == '<');
  istr >> v.x;
  istr >> c;  assert (c == ',');
  istr >> v.y;
  istr >> c;  assert (c == ',');
  istr >> v.z;
  istr >> c;  assert (c == '>');
  return istr;
}

void GiveRainbowColor(double position, unsigned char c[]);
#endif
