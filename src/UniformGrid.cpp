#include "UniformGrid.h"
#include "UniformCell.h"
#include "mesh.h"
#include "hash.h"
#include "triangle.h"
#include "boundingbox.h"
#include "particle.h"
#include "hit.h"
#include <vector>

typedef unsigned int uint;

void UniformGrid::loadMesh( Mesh * mesh , uint d){

  std::cout << "Loading Mesh into uniform grid" << std::endl;

  // get bbox
  const BoundingBox bbox = mesh->getBoundingBox();

  // Settting things inside out box
	double x_range = bbox.getMax().x - bbox.getMin().x;
	double y_range = bbox.getMax().y - bbox.getMin().y;
	double z_range = bbox.getMax().z - bbox.getMin().z;

	dx = x_range / d;
	dy = y_range / d;
	dz = z_range / d;

	division = d;

  // Create the cells where we will store
  for(uint i = d; i < d*d*d; i++)
    uniformCells.push_back(new UniformCell());
  
	// Iterate through all triangles to add them to our mesh
  for ( triangleshashtype::iterator iter = mesh->triangles.begin();
        iter != mesh->triangles.end(); iter++) {
      
      Triangle * t = iter->second;
      insertTriangle( t );

  }

}


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
  // TODO make this function
  // Public function
	assert(false);
}

// Used to get hit object of where a particle hit our mesh
void UniformGrid::collisionDetected( const Particle * p, Hit * h ){
  // TODO make this function
  // Private function
	assert(false);
}

// Get the vector of triangles
std::vector<Triangle * > & UniformGrid::getTriangles(const glm::vec3 & pos){
  return getCell(pos)->read();
}
// Given a point in the space, return a cell
UniformCell * UniformGrid::getCell(const glm::vec3 & position){
  // TODO make this function
	uint i = (int) position.x / dx;
	uint j = (int) position.y / dy;
	uint k = (int) position.z / dz;
  return getCell(i,j,k);
}


// You will you this to insert all triangles into our UG
void UniformGrid::insertTriangle(Triangle * t){

  glm::vec3 a = (*t)[0]->getPos();
  glm::vec3 b = (*t)[1]->getPos();
  glm::vec3 c = (*t)[2]->getPos();    

  // Making a bouding box that holds this triangle
  BoundingBox bb;
  bb.Extend(a);
  bb.Extend(b);
  bb.Extend(c);

  glm::vec3 min = bb.getMin();
  glm::vec3 max = bb.getMax();

  uint i_0 = (uint)(min.x/dx);
  uint j_0 = (uint)(min.y/dy);
  uint k_0 = (uint)(min.z/dz);

  uint i_f = (uint)(max.x/dx);
  uint j_f = (uint)(max.y/dy);
  uint k_f = (uint)(max.z/dz);

  // Create pointer

  // For this range of blocks add in our triangle
  for(uint i =  i_0; i <= i_f; i++)
    for(uint j =  j_0; j <= j_f; j++)
      for(uint k =  k_0; k <= k_f; k++)
        addTriangleToCell(i,j,k,t);
}


// Given the i,j,k return a cell
UniformCell * UniformGrid::getCell(uint i, uint j, uint k){

  // Conversion to 1D array
  int index = i + ( division * j ) + ( division * division * k );
  UniformCell * c = uniformCells[index];
  return c;

}


// Adds a triangle pointer into our grid datastructure
void UniformGrid::addTriangleToCell(uint i , uint j, uint k, Triangle * t){
  UniformCell * c = getCell(i,j,k);
  c->triangles.push_back(t);
}





