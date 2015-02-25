// This is the implementation code for our KDTree Data structure
#include "KDTree.h"

#include <vector>
#include <algorithm>
#include <iostream>

#include "vbo_structs.h"
#include "render_utils.h"
#include "glCanvas.h"

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
void KDTree::update(const std::vector<Particle *> & unsorted, const BoundingBox & b ){

  // Copy boundbox
  bbox = b;

  // Full the heap with empty values
  binary_heap = partPtrVec(unsorted.size(),NULL);

  // If we didnt get any particles
  if(binary_heap.size() == 0)
    return;

  // Create  an unsorted array
  partPtrVec unsorted_copy = unsorted;

  // Create my binary heap rec
  Optimize(0, unsorted_copy.size(), unsorted_copy, 0, 0);
  // std::cout << "Done Building" << std::endl;
  

  for( uint i = 0; i < binary_heap.size(); i++ ){
    if(binary_heap[i] == NULL ){
    
      std::cout << "ERROR INFO HEAP DUMP" << std::endl;
      std::cout << "Size of array A" << unsorted.size() << std::endl;
      for( uint j = 0; j < binary_heap.size(); j++ ){
        std::cout << "binary_heap[" << j <<"] = " << binary_heap[j] << std::endl;
      }
      assert(false);
    
    }
  }

}

void KDTree::Optimize(uint i, uint j, partPtrVec & a, uint8 d, uint hp){
  // Reminders : this will always be changing vector a in place
  // FIXME: Might not work inplace

  std::cout << "Optimize(" << i << ", " << j << ", a, " << d << ", " << hp << ");\n";
  // Do we have a legal hp?p
  if( binary_heap.size() < hp ){
  
    std::cout << "Rejected: " << hp << std::endl;
    return;
  }


  if( i + 1 == j ){ // No need to sort
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
  unsigned int median_index = ((j - i) / 2) + i;
  binary_heap[hp] = a[median_index];

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


// Used to render  for debuging, will render children of  hp
void KDTree::renderKDTree(uint hp, uint8 d, const glm::vec3 & minPt, const glm::vec3 & maxPt) {

  // Setup
  assert( 0 <= hp && hp < binary_heap.size() );
  assert( 0 <= d && d <= 2 );


  glm::vec3 minLeft = minPt;
  glm::vec3 maxLeft;

  glm::vec3 minRight;
  glm::vec3 maxRight = maxPt;

  glm::vec3 point = binary_heap[hp]->getPos();


  if(d == 0){
    // Split along x axis
    
    maxLeft = glm::vec3(point.x, 
                        maxPt.y,
                        maxPt.z);
  
    minRight = glm::vec3(point.x,
                         minPt.y,
                         minPt.z);
  
  }else if( d == 1){
    // Split along y axis
    maxLeft = glm::vec3(maxPt.x, 
                        point.y,
                        maxPt.z);
  
    minRight = glm::vec3(minPt.x,
                         point.y,
                         minPt.z);
  }else{
    // Split along z axis
    maxLeft = glm::vec3(maxPt.x, 
                        maxPt.y,
                        point.z);
  
    minRight = glm::vec3(minPt.x,
                         minPt.y,
                         point.z);
  
  }

  // Render the two boxes created from the split
  // TODO: refactor this so that I do not recreate so many BBOX
  renderBBox(minLeft, maxLeft);
  renderBBox(minRight, maxRight);

  // If you have children left or right recurse on them
  if( hasLeft(hp) )
    renderKDTree(hp*2+1, (d+1)%3, minLeft, maxLeft);

  if( hasRight(hp))
    renderKDTree(hp*2+2, (d+1)%3, minRight, maxRight);

}

void KDTree::renderBBox(const glm::vec3 &A, const glm::vec3 &B) {

  float thickness = 0.0005 * glm::length( B - A );
  glm::vec4 black(0,0,0,1);

  // yz plane
  addEdgeGeometry(tree_verts, tree_tri_indices, glm::vec3(A.x,A.y,A.z), glm::vec3(A.x,A.y,B.z), black,black,thickness,thickness);
  addEdgeGeometry(tree_verts, tree_tri_indices, glm::vec3(A.x,A.y,B.z), glm::vec3(A.x,B.y,B.z), black,black,thickness,thickness);
  addEdgeGeometry(tree_verts, tree_tri_indices, glm::vec3(A.x,B.y,B.z), glm::vec3(A.x,B.y,A.z), black,black,thickness,thickness);
  addEdgeGeometry(tree_verts, tree_tri_indices, glm::vec3(A.x,B.y,A.z), glm::vec3(A.x,A.y,A.z), black,black,thickness,thickness);

  // yz plane forward
  addEdgeGeometry(tree_verts, tree_tri_indices, glm::vec3(B.x,A.y,A.z), glm::vec3(B.x,A.y,B.z), black,black,thickness,thickness);
  addEdgeGeometry(tree_verts, tree_tri_indices, glm::vec3(B.x,A.y,B.z), glm::vec3(B.x,B.y,B.z), black,black,thickness,thickness);
  addEdgeGeometry(tree_verts, tree_tri_indices, glm::vec3(B.x,B.y,B.z), glm::vec3(B.x,B.y,A.z), black,black,thickness,thickness);
  addEdgeGeometry(tree_verts, tree_tri_indices, glm::vec3(B.x,B.y,A.z), glm::vec3(B.x,A.y,A.z), black,black,thickness,thickness);
  
  addEdgeGeometry(tree_verts, tree_tri_indices, glm::vec3(A.x,A.y,A.z), glm::vec3(B.x,A.y,A.z), black,black,thickness,thickness);
  addEdgeGeometry(tree_verts, tree_tri_indices, glm::vec3(A.x,A.y,B.z), glm::vec3(B.x,A.y,B.z), black,black,thickness,thickness);
  addEdgeGeometry(tree_verts, tree_tri_indices, glm::vec3(A.x,B.y,B.z), glm::vec3(B.x,B.y,B.z), black,black,thickness,thickness);
  addEdgeGeometry(tree_verts, tree_tri_indices, glm::vec3(A.x,B.y,A.z), glm::vec3(B.x,B.y,A.z), black,black,thickness,thickness);

}

void KDTree::initializeVBOs(){

  // std::cout << "Before Buffer " << (int)tree_verts_VBO  << std::endl;
  // std::cout << "Buffer " << tree_tri_indices_VBO  << std::endl;
  glGenBuffers(1, &tree_verts_VBO);
  glGenBuffers(1, &tree_tri_indices_VBO);
  // std::cout << "After Buffer " << (int)tree_verts_VBO  << std::endl;
  // std::cout << "Buffer " << (int)tree_tri_indices_VBO  << std::endl;

}

void KDTree::setupVBOs(){

  // std::cout << "Setup Buffer " << (int)tree_verts_VBO  << std::endl;

  tree_verts.clear();
  tree_tri_indices.clear();
  

  if(binary_heap.size() == 0)
    return;

  // setups up all my VBOs
  renderKDTree(0,0,bbox.getMin(), bbox.getMax());

  HandleGLError("enter error verts");
  glBindBuffer(GL_ARRAY_BUFFER,tree_verts_VBO); 
  HandleGLError("leave error verts");
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*tree_verts.size(),&tree_verts[0],GL_STATIC_DRAW); 

  HandleGLError("enter error tri");
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,tree_tri_indices_VBO); 
  HandleGLError("leave error tri");
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(VBOIndexedTri)*tree_tri_indices.size(),&tree_tri_indices[0],GL_STATIC_DRAW);

}

void KDTree::drawVBOs(){

  if(binary_heap.size() == 0)
    return;

  glBindBuffer(GL_ARRAY_BUFFER, tree_verts_VBO);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, tree_tri_indices_VBO);
  glDrawElements(GL_TRIANGLES,
                 tree_tri_indices.size()*3,
                 GL_UNSIGNED_INT, 0);
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);

}

void KDTree::cleanupVBOs(){

  // std::cout << "Someone called cleanupVBOs" << std::endl;
  glDeleteBuffers(1, &tree_verts_VBO);
  glDeleteBuffers(1, &tree_tri_indices_VBO);

}
