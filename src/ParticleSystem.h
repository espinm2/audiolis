#ifndef _PARTICLESYSTEM_H_
#define _PARTICLESYSTEM_H_

/* [x] TODO Add Mesh Ptr to class
 * [x] TODO Expand Constrcutor to take in mesh
 * [x] TODO Create a uniform sphere of dots
 * [x] TODO Alter Load() to take in inputs for simulation
 * [x] TODO Run simulation where spheres project outwards
 * [x] TODO Be able to create that uniform sphere by your cursor
 */

#include <vector>
#include <iostream>
#include <glm/glm.hpp>

#include "glCanvas.h"
#include "boundingbox.h"
#include "vbo_structs.h"


// Forward declaration
class ArgParser;
class Particle;
class BoundingBox;
class Mesh;

class ParticleSystem {


  //TODO move comment over here to header
  public:

    // Constructors
    ParticleSystem(ArgParser *_args, Mesh * _mesh, BoundingBox * _bbox) { 
      args = _args; 
      mesh = _mesh;
      bbox = _bbox; }

    ~ParticleSystem();

    // User interace functions
    void moveCursor(const float & dx, 
        const float & dy, const float & dz );


    // Simulation functions
    void load(); 
    void update();
    bool moveParticle(Particle * p, float timestep);
    void calcMeshCollision(Particle * &p);
    void createInitWave();
    void particleSplit(Particle * &p, std::vector<Particle *> &vec);
    bool shouldSplit(Particle * &p);
    void particleMerge(const Particle * &a, const Particle * &b, Particle * &c);

    // Math function
    double absorbFunc(const std::string & materialName, const double freq);
    
    
    // located in ParticleSystem_render.cpp
    void initializeVBOs();
    void setupVBOs();
    void drawVBOs();
    void cleanupVBOs();


  private:


    // Functions
    void setupParticles();
    void drawParticles();

    void setupCursorPoint();
    void drawCursorPoint();

    // Memebers
    ArgParser * args;
    BoundingBox * bbox;
    Mesh * mesh;
    
    std::vector<Particle *> particles;
    std::vector<Particle *> newParticles;
    glm::vec3 cursor;

    // VBOs Ids
    GLuint particle_verts_VBO;
    GLuint cursor_verts_VBO;

    // Vertices for VBOs
    std::vector<VBOPosNormalColor> particle_verts;
    std::vector<VBOPosNormalColor> cursor_verts;

};

// ===================================================
#endif
