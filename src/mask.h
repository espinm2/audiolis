#ifndef MASK_H
#define MASK_H

#include "particle.h"
#include "render_utils.h"
#include  <vector>


// This is a tempurary datastructure that will create
class Mask{

  public:

    // Constructor
    Mask(Particle * cen, std::vector <Particle *> maskP, int s = 6){

      shape = s;
      maskCenter = cen;
      maskParticles = maskP;

    }

    Mask(){ }
    
    // Standard Accessors
    const std::vector<Particle *> & getMaskParticles() { return maskParticles;}
    Particle * getCenter() { return maskCenter;}
    const std::vector<int> & getCostVector(){ return costVector;}
    int getShape(){ return shape; }

    // Easier to use Accessors
    Particle * getMaskParticle(int i){ return maskParticles[i];}
    const int getCost(int i){ return costVector[i];}
    const int size(){ return size_of_mask; }
    double maskArea(); // TODO

    // Standard setters
    void setCenter( Particle * p){ maskCenter = p; }
    void setMaskParticles( const std::vector<Particle *> pVec ){ maskParticles = pVec;}
    void setCostVector( const std::vector<int> costVec){ costVector = costVec;}
    void setSize( int s ) { size_of_mask = s; }
    void setShape( int s ) {shape = s;}

    // Used to help render this data structure
    void renderOutline( std::vector<VBOPosNormalColor> & outline_verts);
    void renderCost( std::vector<VBOPosNormalColor> & cost_verts);

    // Splits rules for mask
    bool resSpit(std::vector<glm::vec3> & newPartPos);

    // Debug functions
    void debugPrint();
    int getTotalCost();

  private:

    // Pentagon == 5 and Hexagon == 6
    int shape;

    // Center of the mask
    Particle * maskCenter;

    // p1 p2 p3 p4 NULL p6
    // d1 d2 d3 d4 d5 d6

    // Standard convention to use NULL pointer if there is no particle there
    std::vector< Particle *> maskParticles;
    // std::vector< glm::vec3 > delusionalParticleLocations;

    // Standard convention to -1 if there is no cost associated
    std::vector< int > costVector;

    int size_of_mask;

};

#endif // MASK_H
