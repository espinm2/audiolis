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



class UniformGrid {

	public:

	    // Default constructor, requires load to be called
	    UniformGrid(){}

		// Load a mesh into our uniform grid
		void loadMesh( Mesh * mesh, uint division );

		// Used to aid in rendering the mesh for debugging
		void renderGrid( std::vector<VBOPosNormalColor> & buffer);

	    // Prints out the average density  of the mesh to cout
	    void averageDensity();

	    // Get the vector of triangles
	    std::vector<Triangle *> & getTriangles(const glm::vec3 & position);

	    // Rendering functions
	    void initializeVBOs();

	    void setupVBOs();

	    void drawVBOs();

	    void cleanupVBOs();

	private:

		// Given a point in the space, return a cell
		UniformCell * getCell(const glm::vec3 & position);

	    // Given the i,j,k return a cell
	    UniformCell * getCell(uint i, uint j, uint k);

	    // Add Triangle to uniform grid data structure in cell
	    void addTriangleToCell(uint i , uint j, uint k, Triangle * triangle);

		// You will you this to insert all triangles into our UG
		void insertTriangle(Triangle * t);

    	// Reprsentation
		std::vector < UniformCell * > uniformCells; // where our cells are
		double dx,dy,dz; unsigned int division;
		BoundingBox bbox;

		// Bufffers for rendering
	    GLuint uni_verts_VBO;						
	    GLuint uni_tri_indices_VBO;
	    std::vector<VBOPosNormalColor> uni_verts;
	    std::vector<VBOIndexedTri> uni_tri_indices;

};

#endif // UNIFORM_GRID_H
