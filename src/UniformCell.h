#ifndef UNIFORM_CELL_H
#define UNIFORM_CELL_H

#include <vector>

typedef unsigned int uint;

class Mesh;
class UniformCell;
class Triangle;
class Particle;
class Hit;

class UniformCell {

  // Default Constructor
  UniformCell(){}

	public:
		// Where we will store our triangles contained in a cell
		std::vector < Triangle * > triangles;

};

#endif // UNIFORM_CELL_H
