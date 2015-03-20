#include "UniformCell.h"
#include "mesh.h"
#include "hash.h"
#include "triangle.h"
#include "boundingbox.h"
#include "particle.h"
#include "hit.h"
#include <vector>
#define PADDING 0.001

#include "UniformGrid.h"

typedef unsigned int uint;

// Main function called to create our grid
void UniformGrid::loadMesh( Mesh * mesh , uint d){
  std::cout << "Loading Mesh into uniform grid" << std::endl;
  glm::vec3 bmin = mesh->getBoundingBox().getMin();
  glm::vec3 bmax = mesh->getBoundingBox().getMax();
  std::cout << "BBOX_MIN (" << bmin.x << "," << bmin.y << "," << bmin.z<<")\n";
  std::cout << "BBOX_MAX (" << bmax.x << "," << bmax.y << "," << bmax.z<<")\n";

  // Set out bounding box (locked at default)
  glm::vec3 min(-10, -10, -10);
  glm::vec3 max( 10,  10,  10);
  bbox =  BoundingBox(min,max);

  // Setting our dx,dy,dz values
	dx = 20.0 / d;
  dy = 20.0 / d;
  dz = 20.0 / d;

  // Saving how many divisions we have
	division = d;

  // Create the empty cells where we bin our triangles
  for(uint i = 0; i < d*d*d; i++)
    uniformCells.push_back(new UniformCell());
  
	// Iterate through all triangles to add them to our mesh
  for ( triangleshashtype::iterator iter = mesh->triangles.begin();
        iter != mesh->triangles.end(); iter++) {
      insertTriangle( iter->second );
  }
}

// Calculates how many triangles we have located in a bin
void UniformGrid::averageDensity(){

  double sum = 0;
  int n = 0;
  for( UniformCell * cell: uniformCells ){
    if (cell->size() != 0){
      n++; sum += cell->size(); }
  }


  std::cout << "Density of mesh grid object with " 
    << division << " divisions: " << sum / n << std::endl;
}

// Used to aid in rendering the mesh for debugging
void UniformGrid::renderGrid( std::vector<VBOPosNormalColor> & buffer){

  for(uint i = 0; i < division; i++ )
    for(uint j = 0; j < division; j++ )
      for(uint k = 0; k < division; k++ )

}

// Rendering functions 
void UniformGrid::initializeVBOs(){
  glGenBuffers(1, &uni_verts_VBO);
  glGenBuffers(1, &uni_tri_indices_VBO);
}

void UniformGrid::setupVBOs(){

}

void UniformGrid::drawVBOs(){

  glBindBuffer(GL_ARRAY_BUFFER, uni_verts_VBO);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, uni_tri_indices_VBO);
  glDrawElements(GL_TRIANGLES,
                 uni_tri_indices.size()*3,
                 GL_UNSIGNED_INT, 0);
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);

}

void UniformGrid::cleanupVBOs(){

}

// Easy to use function that gives you triangles based on cell
std::vector < Triangle * > & UniformGrid::getTriangles(const glm::vec3 & pos){
  return getCell(pos)->read();
}

// Given a point in the world space, return cell
UniformCell * UniformGrid::getCell(const glm::vec3 & position){

  // Set this point into our space
  glm::vec3 _pos =  position - bbox.getMin();

  // Get where we are in by bin
	uint i = (int) _pos.x / dx;
	uint j = (int) _pos.y / dy;
	uint k = (int) _pos.z / dz;

  return getCell(i,j,k);
}

// Given the i,j,k return a cell
UniformCell * UniformGrid::getCell(uint i, uint j, uint k){
  // Getting the index of our 1D array with conversion
  return uniformCells[i + ( division * j ) + ( division * division * k )];
}

// You will you this to insert all triangles into our UG
void UniformGrid::insertTriangle(Triangle * t){

  // Get points of triangle 
  glm::vec3 a = (*t)[0]->getPos();
  glm::vec3 b = (*t)[1]->getPos();
  glm::vec3 c = (*t)[2]->getPos();

  // Making a bouding box that holds this triangle
  BoundingBox bb; bb.Extend(a); bb.Extend(b); bb.Extend(c);
  glm::vec3 min = bb.getMin(); glm::vec3 max = bb.getMax();

  // Checking to see if the triangle is axis aligned
  glm::vec3 n = t->getNormal();

  /*
  // Adding padding if axis aligned
  glm::vec3 x_axis(1,0,0); glm::vec3 y_axis(0,1,0); glm::vec3 z_axis(0,0,1);
  if( n == x_axis ){min.x -= PADDING; max.x += PADDING; }
  if( n == y_axis ){min.y -= PADDING; max.y += PADDING; }
  if( n == z_axis ){min.z -= PADDING; max.z += PADDING; }
  */

  // Change these from world coordinates to our [-0.5,0.5] World
  min = min - bbox.getMin();
  max = max - bbox.getMin();

  uint i_0 = (int)(min.x/dx);
  uint j_0 = (int)(min.y/dy);
  uint k_0 = (int)(min.z/dz);

  uint i_f = (int)(max.x/dx);
  uint j_f = (int)(max.y/dy);
  uint k_f = (int)(max.z/dz);

  // Create pointer

  // For this range of blocks add in our triangle
  for(uint i =  i_0; i < i_f; i++)
    for(uint j =  j_0; j < j_f; j++)
      for(uint k =  k_0; k < k_f; k++)
        addTriangleToCell(i,j,k,t);
}

// Adds a triangle pointer into our grid datastructure
void UniformGrid::addTriangleToCell(uint i , uint j, uint k, Triangle * t){
  // In range asserts
  assert(0 <= i && i < division);
  assert(0 <= j && j < division);
  assert(0 <= k && k < division);

  getCell(i,j,k)->triangles.push_back(t);
}
