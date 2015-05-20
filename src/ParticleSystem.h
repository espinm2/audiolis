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
#include <fstream>
#include <glm/glm.hpp>

#include "glCanvas.h"
#include "boundingbox.h"
#include "vbo_structs.h"
#include "mask.h"
#include "KDTree.h"
#include "UniformGrid.h"
#include "Stats.h"
#include "sphere.h"

typedef std::vector<std::vector<int>> vMat;
typedef std::vector<Particle *> PartPtrVec;
typedef std::pair<glm::vec3, PartPtrVec> Attractor;
typedef std::vector<Attractor*> AttractorVector;
typedef unsigned int uint;

// Forward declaration
class ArgParser;
class Particle;
class BoundingBox;
class Mesh;
class BVHNode;
class Sphere;

class ParticleSystem {

  friend class Particle; // particle access to mesh

  //TODO move comment over here to header
  public:

    // Constructors
    ParticleSystem(ArgParser *_args, Mesh * _mesh, BoundingBox * _bbox) { 
      args = _args; 
      mesh = _mesh;
      bbox = _bbox; 
    }

    ~ParticleSystem();

    // User interace functions
    void moveCursor(const float & dx, 
        const float & dy, const float & dz );

    // Setter
    void setTimeStep( const double ts){ TIME_STEP = ts; }

    Particle * createParticle(const glm::vec3 & pos, const glm::vec3 & old, 
      const glm::vec3 cen, double watts, double freq, int s);

    // Getter
    unsigned int numParticles(){ return particles.size();}
    double getTimeStep(){ return TIME_STEP;}
    glm::vec3 getCursor(){ return cursor; }

    // Main Simulation Functions
    void load();  // load inital values from args file and meshes
    void update(); // moves, splits, and merges particle in a timestep

    // Events human trigger in simulation
    void createInitWave(); // Creates a sphere of particles
    void stabalizeInitalSphere(); // Reconfigures sound source
    
    // Debug Functions
    void createDebugParticle(); // used for debugging in testing_chamber_1.obj
    void particleSplit(Particle * &p, std::vector<Particle *> &vec);
    void debug(); // Used to test non-visual code
    void collisionDetection(Particle * p);

    // Kuhn–Munkres algorithm based matching
    void munkresMatching(PartPtrVec & partVec, vMat & matchingMat, vMat & costMat,int shape = 6);
    void generateMask(PartPtrVec & conciderForMask, Mask &m);
    void maskFitting(PartPtrVec & conciderForMask, Mask &m, int shape=6);

    // Getting delusional particle locations
    void delusionalParticleLocations(Particle * &cur,
        PartPtrVec &gathered, std::vector<glm::vec3> & output, int shape = 6);

    // Particle Search Functions KD and linear
    void linearGatherParticles(Particle * center, double r, double a, 
        PartPtrVec & result);
    bool linearDuplicateSearch(const glm::vec3 & pos, double th);
    bool linearNewDuplicateSearch(const glm::vec3 & pos, 
        const PartPtrVec & newVec , double th);

    // Miscellaneous <I coulnd't find a group>
    void closeProfiler(){output_profiler_str.close(); }

    // New Update Function Code (tested)
    void generateResSplits(Particle * &cur, PartPtrVec & gathered);
    void mergeSimilarParticles(Particle * &cur, PartPtrVec & gathered);
    void mergeGlobalParticles(double dist);
    void resolveCollisions(Particle * &cur);
    void removeDeadParticles();
    void addNewParticles(); // adds new particles to main vector

    // Relaxation used for inital sphere
    double simulatedannealing(Particle * p, PartPtrVec & gathered); // Moves based on gathered
    void annealing(unsigned int iterations, double prevForce);
    void recompute_collisions();
    glm::vec3 interParticleForce(Particle * & cur, PartPtrVec & partv);

    // New Update Function Code (tested)
    void moveParticle(Particle * & cur); // moves a particle

    // Merge helper functions
    Particle * particlePairMerge(Particle * &a, Particle * &b); 
    Particle * particleVectorMerge(std::vector<Particle *> &vec);

    // Audio based functions
    double absorbFunc(const std::string & materialName, const double freq);

    // Analysis code
    void analyze();
    uint getLowestCostShape(Particle * cur);
    bool getBestFit(Particle * cur, std::vector<glm::vec3> & points);


    // Failed attempt at annealing localized (REMOVE)
    bool maintainDensity(Particle * cur,PartPtrVec & gathered_particles, Attractor * ap);
    bool shouldSplit(Particle * cur,PartPtrVec & gathered_particles);
    double constrainedNudge(Particle * p, PartPtrVec & gathered_particles, glm::vec3 ap_pos);
    void constrainedAnnealing( AttractorVector & av, unsigned int iterations,double prevForce );
    void localAnnealing(unsigned int iterations, double prevForce, 
        std::vector<bool> & fixed, PartPtrVec & gutted_mask_created);
    void prepareMask(PartPtrVec & tmpvec, Mask & m, std::vector<bool> & fixed, PartPtrVec & exposed);

    // This function
    void crunchGaps(const std::vector<double> & sorted_angles, 
        std::vector<double> & gaps);

    // Public render functions used by glCanvas
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

    void setupEdges();
    void drawHappinessVisual();
    
    void setupDelusionalParticles();
    void drawDelusionalParticles();
    void drawDelusionalConnections();

    void setupSphere();
    void drawSphere();

    void setupHalfEdges();
    void drawHalfEdges();

    // Memebers Borrowed from GLCanvas
    ArgParser * args;
    BoundingBox * bbox;
    Mesh * mesh;

    // Memebers unique to this class
    std::vector<Particle *> particles; // Where we store partilces in current iterations
    std::vector<Particle *> newParticles; // Where we put split particles 
    glm::vec3 cursor; // Where the cursor is in world space
    KDTree particle_kdtree; // Where we store particles in a td tre
    BVHNode * root; // We we will store our mesh for fast accesses
    UniformGrid uniform_grid; // Where we store our mesh object fo easy access
    Stats stats; // Where we will keep stats for analysis 
    Sphere sphere; // sphere we used to help us visualize mask

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
    
    unsigned int      ITERATION;
    unsigned int      PARTICLES_PER_M;


    // Gather particles
    double GATHER_DISTANCE; // How far away we gather particles for splits
    double GATHER_ANGLE;    // Directional distance from us  in RAD
    double MERGE_DISTANCE;  // range 0 to this, such as we merge


    unsigned int RELAXATION_MERGE_TRIGGER;

    // Used to save profiling output
    std::ofstream output_profiler_str;
    

    // VBOs Ids
    GLuint particle_verts_VBO;
    GLuint cursor_verts_VBO;
    GLuint velocity_verts_VBO;
    GLuint velocity_tri_indices_VBO;
    GLuint outline_verts_VBO;                                                   // unused 
    GLuint happyness_verts_VBO;                                                 // unused

    GLuint delusional_verts_VBO; 
    GLuint connection_verts_VBO; 

    // Sphere  ids
    GLuint sphere_verts_VBO;
    GLuint sphere_tri_indices_VBO;

    // Half edges
    GLuint half_edges_verts_VBO;

    // Vertices for VBOs
    std::vector<VBOPosNormalColor> particle_verts;
    std::vector<VBOPosNormalColor> cursor_verts;
    std::vector<VBOPosNormalColor> velocity_verts;
    std::vector<VBOIndexedTri> velocity_tri_indices;
    std::vector<VBOPosNormalColor> outline_verts;
    std::vector<VBOPosNormalColor> happyness_verts;

    // Experimental
    std::vector<VBOPosNormalColor> delusional_verts; // TODO
    std::vector<VBOPosNormalColor> connection_verts; // TODO

    // Sphere
    std::vector<VBOPosNormalColor> sphere_verts;
    std::vector<VBOIndexedTri>     sphere_tri_indices;

    // Half edges
    std::vector<VBOPosNormalColor> half_edges_verts;

};



// Utility class to sort particles by their distance relative to another
class particleCMP{

  glm::vec3 c; // center particle we will use to compare

public:

  particleCMP(const glm::vec3 & center){
    c = center;
  }

  bool operator() (Particle * lhs ,Particle * rhs) const{
    glm::vec3 l = lhs->getOldPos(); glm::vec3 r = rhs->getOldPos();
    double dist_l_squared = pow(c.x - l.x ,2) + pow(c.y - l.y ,2) + pow(c.z - l.z ,2);
    double dist_r_squared = pow(c.x - r.x ,2) + pow(c.y - r.y ,2) + pow(c.z - r.z ,2);

    return dist_l_squared < dist_r_squared;
  }
};

#endif
