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
#include "vbo_structs.h"
#include <GLFW/glfw3.h>

class ArgParser;
class Particle;

class ParticleSystem {


  public:

    // Constructors
    ParticleSystem(ArgParser *_args) { args = _args; }
    ~ParticleSystem();

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

    // Memebers
    ArgParser * args;
    std::vector<Particle *> particles;

    // VBOs
    GLuint particle_verts_VBO;
    std::vector<VBOPosNormalColor> particle_verts;
    

};

// ===================================================
#endif
