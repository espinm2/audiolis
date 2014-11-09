#ifndef _PARTICLESYSTEM_H_
#define _PARTICLESYSTEM_H_


/* [x] TODO Add Mesh Ptr to class
 * [x] TODO Expand Constrcutor to take in mesh
 * [ ] TODO Create a uniform sphere of dots
 * [ ] TODO Alter Load() to take in inputs for simulation
 * [ ] TODO Run simulation where spheres project outwards
 * [ ] TODO Be able to create that uniform sphere by your cursor
 * [ ] TODO Work on Collision Detection, Make triangles turn red if hit
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
    void moveCursor(const double & dx, 
        const double & dy, const double & dz );


    // Simulation functions
    void load(); 
    void update();
    
    // located in ParticleSystem_render.cpp
    void initializeVBOs();
    void setupVBOs();
    void drawVBOs();
    void cleanupVBOs();


  private:

    // Fuctions for Faking Depth
    int getGLPointSize(const Vec3f & point);

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
