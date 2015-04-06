#ifndef __ARG_PARSER_H__
#define __ARG_PARSER_H__

#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <string>
#include "MersenneTwister.h"

// ================================================================================
// ================================================================================

inline void separatePathAndFile(const std::string &input, std::string &path, std::string &file) {
  // we need to separate the filename from the path
  // (we assume the vertex & fragment shaders are in the same directory)
  // first, locate the last '/' in the filename
  size_t last = std::string::npos;  
  while (1) {
    int next = input.find('/',last+1);
    if (next != std::string::npos) { 
      last = next;
      continue;
    }
    next = input.find('\\',last+1);
    if (next != std::string::npos) {
      last = next;
      continue;
    }
    break;
  }
  if (last == std::string::npos) {
    // if there is no directory in the filename
    file = input;
    path = ".";
  } else {
    // separate filename & path
    file = input.substr(last+1,input.size()-last-1);
    path = input.substr(0,last);
  }
}

// ======================================================================
// Class to collect all the high-level rendering parameters controlled
// by the command line or the keyboard input
// ======================================================================


class ArgParser {

public:

  ArgParser() { DefaultValues(); }

  ArgParser(int argc, char *argv[]) {
    DefaultValues();
    // parse the command line arguments
    for (int i = 1; i < argc; i++) {
      if (std::string(argv[i]) == std::string("-input") || 
          std::string(argv[i]) == std::string("-i")) {
        i++; assert (i < argc); 
        separatePathAndFile(argv[i],path,input_file);

      } else if (!strcmp(argv[i],"-shader")) {
        shader_filename = std::string(argv[i]);

      } else if (!strcmp(argv[i],"-size")) {
        i++; assert (i < argc); 
        width = height = atoi(argv[i]);

      } else if (!strcmp(argv[i],"-d")) {
        i++; assert (i < argc); 
        division = atoi(argv[i]);

      } else if (!strcmp(argv[i],"-init_particles")) {
        i++; assert (i < argc); 
        num_init_particles = atoi(argv[i]);

      } else if (!strcmp(argv[i],"-timestep")) {
        i++; assert (i < argc); 
        timestep = atof(argv[i]);

      } else if (!strcmp(argv[i],"-fps")) {
        i++; assert (i < argc); 
        fps_cap = atoi(argv[i]);

      } else if (!strcmp(argv[i],"-profile")) {
        profile = true;

      } else {
        printf ("whoops error with command line argument %d: '%s'\n",i,argv[i]);
        assert(0);
      }
    }
  }

  void DefaultValues() {

    input_file = "";
    path = "";
    shader_filename = "../src/shader";

    width = 800;
    height = 600;

    geometry = false;
    bounding_box = false;
    gouraud_normals = false;

    timer = 0.0;
    timestep = -1;
    num_init_particles = 1000;
    
    animate = false;
    wireframe = false;
    render_top = true;
    source_type = 0;
    viz_type = 3;
    direction = false;

    wall_material = 0;
    floor_material = 0;
    ceiling_material = 0;
    absorber = false;

    printcusorpos = true;
    render_outline = false;
    render_edges = false;
    render_mask = 0;

    output_file = "merge_profiling_data.txt";
    profile = false;
    setupInitParticles = true;
    kdtree_render = false;
    ugrid_render = false;

    division = 10;
    fps_cap = 30;


  }

  // ==============
  // REPRESENTATION
  // all public! (no accessors)

  std::string input_file;
  std::string path;
  std::string shader_filename;

  std::string output_file;

  int width;
  int height;

  bool geometry;        
  bool bounding_box;
  bool gouraud_normals;

  float timer; // acts like global var for FPS counter
  double timestep; // passed to ParticleSystem.cpp
  unsigned int num_init_particles; // passed to ParticleSystem.cpp

  // Toggles visualizations also acts like global varibles
  bool animate;       
  bool wireframe;
  bool render_top;
  int source_type;
  int viz_type;
  bool direction;
  bool printcusorpos;

  int wall_material; // toggle what the walls are made of
  int floor_material; //  Toggle what the floor is made of
  int ceiling_material; // Toggle what the ceilings are made of
  bool absorber;  // Toggle if there is an absorber in the room

  bool render_outline;
  bool render_edges;

  int render_mask;

  bool setupInitParticles;
  bool profile;

  bool kdtree_render;
  bool ugrid_render;

  int division;
  int fps_cap;

  // Material Key
  // wall_material 0 => bricked_wall; 
  //               1 => concrete_wall; 
  //               2 => ceramic_wall;
  //
  // floor_material 0 =>  pvc_floor;
  //                1 => carpated_floor;
  //
  // ceiling set to plaster_ceiling
  // glass is set to double plated glass
  

  // Random number generator
  MTRand randomGen;
};


#endif
