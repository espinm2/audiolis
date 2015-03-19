#ifndef UNIFORM_GRID_H
#define UNIFORM_GRID_H

#include "vbo_structs.h"

/*
Data structure will be used to  reduce the amount of time required to
handle collisions between particles and the walls our space.
This is a static step that will only happen upon creation.
*/

#include <vector>

typedef unsigned int uint;


// Increases compile time
class Mesh;
class UniformCell;
class BoundingBox;
class Triangle;
class Particle;
class Hit;



class UniformGrid {

	public:

		// Given a mesh will break apart the scene into grids
		UniformGrid(uint division, const BoundingBox * bbox);

		// Used to aid in rendering the mesh for debugging
		void renderGrid( std::vector<VBOPosNormalColor> & buffer);

		// Used to get hit object of where a particle hit our mesh
		void collisionDetected( const Particle * p, Hit * h );


	private:

		// Load a mesh into our uniform grid
		void loadMesh( Mesh * mesh );

		// Given a point in the space, return a cell
		UniformCell * getCell(const glm::vec3 & position);

    // Given the i,j,k return a cell
    UniformCell * getCell(uint i, uint j, uint k);

    // Add Triangle to uniform grid data structure in cell
    void addTriangleToCell(uint i , uint j, uint k, Triangle * triangle);

		// You will you this to insert all triangles into our UG
		void insertTriangle(Triangle * t);
    // Reprsentation
    
		// Where are store our UniformCells;
		std::vector < UniformCell * > uniformCells;

		// Size our cells
		double dx,dy,dz;
		uint division;

};

#endif // UNIFORM_GRID_H
