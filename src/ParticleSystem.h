#ifndef _PARTICLESYSTEM_H_
#define _PARTICLESYSTEM_H_


/* [x] TODO Add Mesh Ptr to class
 * [x] TODO Expand Constrcutor to take in mesh
 * [x] TODO Create a uniform sphere of dots
 * [x] TODO Alter Load() to take in inputs for simulation
 * [x] TODO Run simulation where spheres project outwards
 * [x] TODO Be able to create that uniform sphere by your cursor
 * [ ] TODO Work on Collision Detection, Make triangles turn red if hit
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
    void moveParticle(Particle * &p);
    
    
    // located in ParticleSystem_render.cpp
    void initializeVBOs();
    void setupVBOs();
    void drawVBOs();
    void cleanupVBOs();


  private:

    // Fuctions for Faking Depth
    int getGLPointSize(const glm::vec3 & point);

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
