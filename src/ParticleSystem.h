#ifndef _PARTICLESYSTEM_H_
#define _PARTICLESYSTEM_H_


/*
 *
 * [ ] TODO Add Mesh Ptr to class
 * [ ] TODO Expand Constrcutor to take in mesh
 * [ ] TODO Alter Load() to take in inputs for simulation
 *
 */
#include <vector>
#include <iostream>
#include "glCanvas.h"
#include "boundingbox.h"

#include "vbo_structs.h"
#include "vectors.h"



// Forward declaration
class ArgParser;
class Particle;
class BoundingBox;
class Mesh;

class ParticleSystem {


  public:

    // Constructors
    ParticleSystem(ArgParser *_args, Mesh * _mesh, BoundingBox * _bbox) { 
      args = _args; 
      mesh = _mesh; 
      bbox = _bbox; }

    ~ParticleSystem();

    // User interace functions
    void moveCursor(const double & dx, const double & dy, const double & dz );


    // Simulation functions
    void load(); 
    void update();
    
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
    Vec3f cursor;

    // VBOs Ids
    GLuint particle_verts_VBO;
    GLuint cursor_verts_VBO;

    // Vertices for VBOs
    std::vector<VBOPosNormalColor> particle_verts;
    std::vector<VBOPosNormalColor> cursor_verts;
    

};

// ===================================================
#endif
