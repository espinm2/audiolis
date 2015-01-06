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

    // Setter
    void setTimeStep( const double ts){ TIME_STEP = ts; }

    // Simulation functions
    void load(); 
    void debug();
    void update();
    bool moveParticle(Particle * p, double timestep);
    void calcMeshCollision(Particle * &p);
    void createInitWave();
    void particleSplit(Particle * &p, std::vector<Particle *> &vec);
    bool shouldSplit(Particle * &p);

    // Two merge particles
    Particle * particlePairMerge(Particle * &a, Particle * &b);
    Particle * particleVectorMerge(std::vector<Particle *> &vec);

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

    void setupVelocityVisual();
    void drawVelocityVisual();

    // Memebers
    ArgParser * args;
    BoundingBox * bbox;
    Mesh * mesh;


    // Simuation Important Varibles
    double            TIME_STEP;  // how much time is passed in seconds
    float             VELOCITY_OF_MEDIUM; // velocity of air in m/s

    double            RADIUS_INIT_SPHERE; // radius of source sphere
    unsigned int      NUM_INIT_PARTICLES; // number of initial particles
    double            MIN_WATTAGE; // when a particle falls below power threshold
    unsigned int      MAX_ITERATIONS; // how many iterations before you kill particles

    // For splits, could be removed
    double            RADIUS_PARTICLE_WAVE; // radius of cluster created upon splits
    unsigned int      SPLIT_AMOUNT; // how many created upon splits, excluding center
    
    std::vector<Particle *> particles; // Where we store partilces in current iterations
    std::vector<Particle *> newParticles; // Where we put split particles 

    glm::vec3 cursor; // Where the cursor is in world space

    // VBOs Ids
    GLuint particle_verts_VBO;
    GLuint cursor_verts_VBO;
    GLuint velocity_verts_VBO;
    GLuint velocity_tri_indices_VBO;

    // Vertices for VBOs
    std::vector<VBOPosNormalColor> particle_verts;
    std::vector<VBOPosNormalColor> cursor_verts;
    std::vector<VBOPosNormalColor> velocity_verts;
    std::vector<VBOIndexedTri> velocity_tri_indices;

};

// ===================================================
#endif
