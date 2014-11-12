#ifndef __ARG_PARSER_H__
#define __ARG_PARSER_H__

#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <string>

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
    geometry = true;
    bounding_box = false;
    gouraud_normals = false;
    timer = 0.0;
    animate = false;
    wireframe = false;
    render_top = false;
    timestep = 0.001;
  }

  // ==============
  // REPRESENTATION
  // all public! (no accessors)

  std::string input_file;
  std::string path;
  std::string shader_filename;
  int width;
  int height;
  bool geometry;
  bool bounding_box;
  bool gouraud_normals;
  float timer;
  bool animate;
  bool wireframe;
  bool render_top;
  float timestep;
};


#endif
