#ifndef MASK_H
#define MASK_H

#include "particle.h"
#include "render_utils.h"
#include  <vector>


// This is a tempurary datastructure that will create
class Mask{

  public:


    // Constructor
    Mask(Particle * cen, std::vector <Particle *> maskP){
      maskCenter = cen;
      maskParticles = maskP; }

    Mask(){ }
    
    // Standard Accessors
    const std::vector<Particle *> & getMaskParticles() { return maskParticles;}
    const Particle * getCenter() { return maskCenter;}
    const std::vector<int> & getCostVector(){ return costVector;}

    // Easier to use Accessors
    const Particle * getMaskParticle(int i){ return maskParticles[i];}
    const int getCost(int i){ return costVector[i];}
    const int size(){ return size_of_mask; }

    // Standard setters
    void setCenter( Particle * p){ maskCenter = p; }
    void setMaskParticles( const std::vector<Particle *> pVec ){ maskParticles = pVec;}
    void setCostVector( const std::vector<int> costVec){ costVector = costVec;}
    void setSize( int s ) { size_of_mask = s; }

    // Used to help render this data structure
    void renderOutline( std::vector<VBOPosNormalColor> & outline_verts);
    void renderCost( std::vector<VBOPosNormalColor> & cost_verts);

    // Splits rules for mask
    bool resSpit(std::vector<glm::vec3> & newPartPos);
    void debugPrint();

  private:

    // Center of the mask
    Particle * maskCenter;

    // Standard convention to use NULL pointer if there is no particle there
    std::vector< Particle *> maskParticles;

    // Standard convention to -1 if there is no cost associated
    std::vector< int > costVector;

    int size_of_mask;

};

#endif // MASK_H
