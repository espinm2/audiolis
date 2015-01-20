// Graphics Library Includes
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// for sleep 	
#if defined(_WIN32)
#include <windows.h>
#else
#include <unistd.h>
#endif

#include "argparser.h"
#include "glCanvas.h"
#include "camera.h"
#include <time.h>

// ====================================================================

int main(int argc, char *argv[]) {

  // parse the command line arguments
  ArgParser args(argc, argv);
  GLCanvas::initialize(&args); 


  // Creates background color behind all objects (SKY)
  // glClearColor(0.8,0.9,1.0,0.0);
  // glClearColor(176.0/256,196.0/256,222.0/256,0.0);
  glClearColor(0.3,0.3,0.3,1.0);
  
  // If enabled do depth comparisons and update the depth buffer
  glEnable(GL_DEPTH_TEST);

  // Specify the value used for depth buffer comparisons
  glDepthFunc(GL_LESS); 
  glDisable(GL_CULL_FACE);

  // Single Instance Calculations /////
  
  // Bounding box
  glm::vec3 center;
  GLCanvas::bbox.getCenter(center);

  // Scaling factor
  double maxDim = GLCanvas::bbox.maxDim();
  float scaleFactor = 2.0 / float(maxDim);
  glm::mat4 myScalingMatrix = glm::scale(glm::vec3(scaleFactor,scaleFactor,scaleFactor));
  
  // While we don't close this window
  while (!glfwWindowShouldClose(GLCanvas::window))  {
    
    // Clear the buffers, color info + depth info
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Run this program GLCanvas
    glUseProgram(GLCanvas::programID);

    // Create camera (VIEW TRANSFORMATIONS)
    GLCanvas::camera->glPlaceCamera();
    glm::mat4 myTranslateMatrix = glm::translate(-center);
    glm::mat4 ModelMatrix = myScalingMatrix*myTranslateMatrix;
    
    // Build the matrix to position the camera based on keyboard and mouse input
    glm::mat4 ProjectionMatrix = GLCanvas::camera->getProjectionMatrix();
    glm::mat4 ViewMatrix = GLCanvas::camera->getViewMatrix();

    // Callings the draw
    GLCanvas::drawVBOs(ProjectionMatrix,ViewMatrix,ModelMatrix);

    // Calling the simulations
    GLCanvas::animate();

    // Swap buffers
    glfwSwapBuffers(GLCanvas::window);
    fflush(stdout);
    glfwPollEvents();  
    fflush(stdout);

#if defined(_WIN32)
  Sleep(10);
#else
  usleep(10);
#endif

  }
  
  GLCanvas::cleanupVBOs();
  glDeleteProgram(GLCanvas::programID);
  
  // Close OpenGL window and terminate GLFW
  glfwDestroyWindow(GLCanvas::window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

// ====================================================================
// ====================================================================
