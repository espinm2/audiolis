// This is the implementation code for our KDTree Data structure

#include "KDTree.h"
#include "particle.h"
#include "vbo_structs.h"
#include <vector>
#include <algorithm>

typedef unsigned int uint;
typedef short unsigned int uint8;
typedef std::vector<Particle *>  partPtrVec;

// Nonmemeber comparsion functions used by Optimize
bool compare_x_pos(Particle * a, Particle * b){
  return a->getPos().x < b->getPos().x; }

bool compare_y_pos(Particle * a, Particle * b){
  return a->getPos().y < b->getPos().y; }

bool compare_z_pos(Particle * a, Particle * b){
  return a->getPos().z < b->getPos().z; }


// Constructor where you have all the elments of the KD tree
KDTree::KDTree(const std::vector<Particle *> & unsorted ){

  // Full the heap with empty values
  binary_heap = partPtrVec(unsorted.size(),NULL);

  // Create  an unsorted array
  partPtrVec unsorted_copy = unsorted;

  // Create my binary heap rec
  Optimize(0, unsorted_copy.size(), unsorted_copy, 0, 0);

}

void KDTree::Optimize(uint i, uint j, partPtrVec & a, uint8 d, uint hp){
  // Reminders : this will always be changing vector a in place
  // FIXME: Might not work inplace

  // Make sure my heap root is correct
  assert(0 <= hp && hp < binary_heap.size());
  assert(i <= j); 

  if(i == j){ // No need to sort
    binary_heap[hp] = a[i]; // cpy pointer for heap
    return;
  }

  // Sort my list depending on my descriminator
  switch(d){

    case 0:
      std::sort(a.begin() + i, a.begin()+j, compare_x_pos);
      break;
    case 1:
      std::sort(a.begin() + i, a.begin()+j, compare_y_pos);
      break;
    case 2:
      std::sort(a.begin() + i, a.begin()+j, compare_z_pos);
      break;
    default:
      assert(false);
  }
  
  // find median
  unsigned int median_index = (j - i) / 2;
  binary_heap[hp] = a[median_index]; // cpy pointer for heap

  d = (d + 1)  % 3; // update descriminator

  // Left child & Right child recur
  Optimize(i             , median_index, a, d, hp * 2 + 1 );
  Optimize(median_index+1, j           , a, d, hp * 2 + 2 );

}
// Used to get the children
const Particle *  KDTree::rightChild(unsigned int index) const{
  return (binary_heap[index*2 + 2]);
}

const Particle *  KDTree::leftChild(unsigned int index)  const{
  return (binary_heap[index*2 + 1]);
}

void KDTree::rightChild(uint index, Particle * p){
  uint i = index * 2 + 2;
  assert(0 <= index && i < binary_heap.size()); // sanity check
  binary_heap[i] = p; // copied the pointer p into my heap

}

void KDTree::leftChild(uint index, Particle * p){
  uint i = index * 2 + 1;
  assert(0 <= index && i < binary_heap.size()); // sanity check
  binary_heap[i] = p; // copied the pointer p into my heap
}


// Used to render  for debuging
void KDTree::renderKDTree( std::vector<VBOPosNormalColor> & outline_verts) const{


}


