#include "glCanvas.h"
#include "argparser.h"
#include "camera.h"

#include <unistd.h>

#if _WIN32
#include <Winsock2.h>
#else
#include <sys/time.h>
#endif

#include "mesh.h"
#include "utils.h"

// ========================================================
// static variables of GLCanvas class

ArgParser* GLCanvas::args = NULL;
Camera* GLCanvas::camera = NULL;

Mesh* GLCanvas::mesh = NULL;
BoundingBox GLCanvas::bbox;
GLFWwindow* GLCanvas::window = NULL;

// mouse position
int GLCanvas::mouseX = 0;
int GLCanvas::mouseY = 0;
// which mouse button
bool GLCanvas::leftMousePressed = false;
bool GLCanvas::middleMousePressed = false;
bool GLCanvas::rightMousePressed = false;
// current state of modifier keys
bool GLCanvas::shiftKeyPressed = false;
bool GLCanvas::controlKeyPressed = false;
bool GLCanvas::altKeyPressed = false;
bool GLCanvas::superKeyPressed = false;

GLuint GLCanvas::render_VAO;

GLuint GLCanvas::ViewMatrixID;
GLuint GLCanvas::ModelMatrixID;
GLuint GLCanvas::LightID;
GLuint GLCanvas::MatrixID;
GLuint GLCanvas::programID;
GLuint GLCanvas::whichshaderID;
GLuint GLCanvas::colormodeID;


// ========================================================
// Initialize all appropriate OpenGL variables, set
// callback functions, and start the main event loop.
// This function will not return but can be terminated
// by calling 'exit(0)'
// ========================================================

void GLCanvas::initialize(ArgParser *_args) {

  args = _args;
  mesh = new Mesh(args);
  mesh->Load();
  bbox.Set(mesh->getBoundingBox());

  glfwSetErrorCallback(error_callback);

  // Initialize GLFW
  if( !glfwInit() ) {
    std::cerr << "ERROR: Failed to initialize GLFW" << std::endl;
    exit(1);
  }
  
  // We will ask it to specifically open an OpenGL 3.2 context
  glfwWindowHint(GLFW_SAMPLES, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  // Create a GLFW window
  window = glfwCreateWindow(args->width,args->height, "OpenGL viewer", NULL, NULL);
  if (!window) {
    std::cerr << "ERROR: Failed to open GLFW window" << std::endl;
    glfwTerminate();
    exit(1);
  }
  glfwMakeContextCurrent(window);
  HandleGLError("in glcanvas first");

  // Initialize GLEW
  glewExperimental = true; // Needed for core profile
  if (glewInit() != GLEW_OK) {
    std::cerr << "ERROR: Failed to initialize GLEW" << std::endl;
    glfwTerminate();
    exit(1);
  }

  // there seems to be a "GL_INVALID_ENUM" error in glewInit that is a
  // know issue, but can safely be ignored
  HandleGLError("after glewInit()",true);

  std::cout << "-------------------------------------------------------" << std::endl;
  std::cout << "OpenGL Version: " << (char*)glGetString(GL_VERSION) << '\n';
  std::cout << "-------------------------------------------------------" << std::endl;

  // Initialize callback functions
  glfwSetCursorPosCallback(GLCanvas::window,GLCanvas::mousemotionCB);
  glfwSetMouseButtonCallback(GLCanvas::window,GLCanvas::mousebuttonCB);
  glfwSetKeyCallback(GLCanvas::window,GLCanvas::keyboardCB);

  programID = LoadShaders(  args->path+"/"+args->shader_filename+".vs",
                            args->path+"/"+args->shader_filename+".fs");

  GLCanvas::initializeVBOs();
  GLCanvas::setupVBOs();

  // ===========================
  // initial placement of camera 
  // look at an object scaled & positioned to just fit in the box (-1,-1,-1)->(1,1,1)
  glm::vec3 camera_position = glm::vec3(1,3,8);
  glm::vec3 point_of_interest = glm::vec3(0,0,0);
  glm::vec3 up = glm::vec3(0,1,0);
  float angle = 20.0;
  camera = new PerspectiveCamera(camera_position, point_of_interest, up, angle);
  camera->glPlaceCamera(); 

  HandleGLError("finished glcanvas initialize");
}




void GLCanvas::animate(){
#if _WIN32    
  static int last_tick_count;
  static int last_fps_count;
  static int frames = -1;
  if (frames == -1) {    
    last_fps_count = last_tick_count = GetTickCount();
    frames = 0;
  }
  if (args->animate) {
    frames++;
    int this_tick_count = GetTickCount();
    // difference in milliseconds
    double diff = 0.001*(this_tick_count-last_tick_count);
    last_tick_count = this_tick_count;
    args->timer += diff;
    double diff_fps_time = 0.001*(this_tick_count-last_fps_count);
    if (diff_fps_time > 1.00) {      
      float fps = frames / float(diff_fps_time);
      std::cout << "fps: " << fps << std::endl;
      frames = 0;
      last_fps_count = this_tick_count;
    }
    mesh->setupVBOs();
  } else {
    last_tick_count = last_fps_count = GetTickCount();
  }
#else
  static timeval last_time;
  static timeval last_fps_time;
  static int frames = -1;
  if (frames == -1) {
    gettimeofday(&last_time,NULL);
    last_fps_time = last_time;
    frames = 0;
  }
  if (args->animate) {
    frames++;
    timeval this_time;
    gettimeofday(&this_time,NULL);
    // compute the difference from last time
    double diff = this_time.tv_sec - last_time.tv_sec + 
      0.000001 * (this_time.tv_usec - last_time.tv_usec);
    double diff_fps_time = this_time.tv_sec - last_fps_time.tv_sec + 
      0.000001 * (this_time.tv_usec - last_fps_time.tv_usec);
    last_time = this_time;
    // print out stats on the FPS occasionally
    if (diff_fps_time > 1.00) {      
      float fps = frames / float(diff_fps_time);
      std::cout << "fps: " << fps << std::endl;
      frames = 0;
      last_fps_time = this_time;
    }
    args->timer += diff;
    mesh->setupVBOs();
  } else {
    gettimeofday(&last_time, NULL);
    last_fps_time = last_time;
    frames = 0;
    // if we aren't animating the light source, avoid busy-waiting!
    //usleep (100);
  }
#endif
  
}


void GLCanvas::initializeVBOs(){
  HandleGLError("enter initilizeVBOs()");
  glGenVertexArrays(1, &render_VAO);
  glBindVertexArray(render_VAO);
  GLCanvas::MatrixID = glGetUniformLocation(GLCanvas::programID, "MVP");
  GLCanvas::LightID = glGetUniformLocation(GLCanvas::programID, "LightPosition_worldspace");
  GLCanvas::ViewMatrixID = glGetUniformLocation(GLCanvas::programID, "V");
  GLCanvas::ModelMatrixID = glGetUniformLocation(GLCanvas::programID, "M");
  GLCanvas::colormodeID = glGetUniformLocation(GLCanvas::programID, "colormode");
  GLCanvas::whichshaderID = glGetUniformLocation(GLCanvas::programID, "whichshader");
  mesh->initializeVBOs();
  HandleGLError("leaving initilizeVBOs()");
}


void GLCanvas::setupVBOs(){
  HandleGLError("enter GLCanvas::setupVBOs()");
  assert (mesh != NULL);
  mesh->setupVBOs();
  HandleGLError("leaving GLCanvas::setupVBOs()");
}


void GLCanvas::drawVBOs(const glm::mat4 &ProjectionMatrix,const glm::mat4 &ViewMatrix,const glm::mat4 &ModelMatrix){
  HandleGLError("enter GlCanvas::drawVBOs()");

  // prepare data to send to the shaders
  glm::mat4 MVP = ProjectionMatrix * ViewMatrix * ModelMatrix;
  glm::vec3 lightPos = mesh->LightPosition();
  glm::vec4 lightPos2 = glm::vec4(lightPos.x,lightPos.y,lightPos.z,1);
  lightPos2 = ModelMatrix * lightPos2;
  glUniform3f(GLCanvas::LightID, lightPos2.x, lightPos2.y, lightPos2.z);

  glUniformMatrix4fv(GLCanvas::MatrixID, 1, GL_FALSE, &MVP[0][0]);
  glUniformMatrix4fv(GLCanvas::ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
  glUniformMatrix4fv(GLCanvas::ViewMatrixID, 1, GL_FALSE, &ViewMatrix[0][0]);

  mesh->drawVBOs();
  HandleGLError("leaving GlCanvas::drawVBOs()");
}


void GLCanvas::cleanupVBOs(){
  bbox.cleanupVBOs();
  mesh->cleanupVBOs();
}


// ========================================================
// Callback function for mouse click or release
// ========================================================

void GLCanvas::mousebuttonCB(GLFWwindow *window, int which_button, int action, int mods) {
  // store the current state of the mouse buttons
  if (which_button == GLFW_MOUSE_BUTTON_1) {
    if (action == GLFW_PRESS) {
      leftMousePressed = true;
    } else {
      assert (action == GLFW_RELEASE);
      leftMousePressed = false;
    }
  } else if (which_button == GLFW_MOUSE_BUTTON_2) {
    if (action == GLFW_PRESS) {
      rightMousePressed = true;
    } else {
      assert (action == GLFW_RELEASE);
      rightMousePressed = false;
    }
  } else if (which_button == GLFW_MOUSE_BUTTON_3) {
    if (action == GLFW_PRESS) {
      middleMousePressed = true;
    } else {
      assert (action == GLFW_RELEASE);
      middleMousePressed = false;
    }
  }
}	

// ========================================================
// Callback function for mouse drag
// ========================================================

void GLCanvas::mousemotionCB(GLFWwindow *window, double x, double y) {

  // camera controls that work well for a 3 button mouse
  if (!shiftKeyPressed && !controlKeyPressed && !altKeyPressed) {
    if (leftMousePressed) {
      camera->rotateCamera(mouseX-x,mouseY-y);
    } else if (middleMousePressed)  {
      camera->truckCamera(mouseX-x, y-mouseY);
    } else if (rightMousePressed) {
      camera->dollyCamera(mouseY-y);
    }
  }

  if (leftMousePressed || middleMousePressed || rightMousePressed) {
    if (shiftKeyPressed) {
      camera->zoomCamera(mouseY-y);
    }
    // allow reasonable control for a non-3 button mouse
    if (controlKeyPressed) {
      camera->truckCamera(mouseX-x, y-mouseY);    
    }
    if (altKeyPressed) {
      camera->dollyCamera(y-mouseY);    
    }
  }
  mouseX = x;
  mouseY = y;
}

// ========================================================
// Callback function for keyboard events
// ========================================================

void GLCanvas::keyboardCB(GLFWwindow* window, int key, int scancode, int action, int mods) {
  // store the modifier keys
  shiftKeyPressed = (GLFW_MOD_SHIFT & mods);
  controlKeyPressed = (GLFW_MOD_CONTROL & mods);
  altKeyPressed = (GLFW_MOD_ALT & mods);
  superKeyPressed = (GLFW_MOD_SUPER & mods);
  // non modifier key actions

  if (key == GLFW_KEY_ESCAPE || key == 'q' || key == 'Q') {
    glfwSetWindowShouldClose(GLCanvas::window, GL_TRUE);
  }

  // other normal ascii keys...
  if ( (action == GLFW_PRESS || action == GLFW_REPEAT) && key < 256) {
    switch (key) {
    case 'b': case 'B':
      args->bounding_box = !args->bounding_box;
      mesh->setupVBOs();
      break;
    case 'g': case 'G':
      args->geometry = !args->geometry;
      mesh->setupVBOs();
      break;
    case 'm': case 'M':
      args->mirror = !args->mirror;
      mesh->setupVBOs();
      break;
    case 'r': case 'R':
      args->reflected_geometry = !args->reflected_geometry;
      mesh->setupVBOs();
      break;
    case 's': case 'S':
      args->shadow = !args->shadow;
      mesh->setupVBOs();
      break;
    case 'e': case 'E':
      args->silhouette_edges = !args->silhouette_edges;
      mesh->setupVBOs();
      break;
    case 'p': case 'P':
      args->shadow_polygons = !args->shadow_polygons;
      mesh->setupVBOs();
      break;
    case 'a': case 'A':
      args->animate = !args->animate;
      break;
    case ' ':
      args->timer += 0.1;
      mesh->setupVBOs();
      break;
    case 'n': case 'N':
      args->gouraud_normals = !args->gouraud_normals;
      mesh->setupVBOs();
      break;
    case 'x': case 'X':
      args->whichshader = (args->whichshader+1)%4;
      std::cout << "CURRENT SHADER: " << args->whichshader;
      if (args->whichshader == 0) { std::cout << " none" << std::endl; }
      if (args->whichshader == 1) { std::cout << " checkerboard" << std::endl; }
      if (args->whichshader == 2) { std::cout << " orange" << std::endl; }
      if (args->whichshader == 3) { std::cout << " wood" << std::endl; }
      mesh->setupVBOs();
      break;
    case 'l' : case 'L':
      //LoadCompileLinkShaders();
      mesh->setupVBOs();
      break;
    case 'q':  case 'Q':
      // quit
      glfwSetWindowShouldClose(GLCanvas::window, GL_TRUE);
      break;
    default:
      std::cout << "UNKNOWN KEYBOARD INPUT  '" << (char)key << "'" << std::endl;
    }
    setupVBOs();
  }
}


// ========================================================
// Load the vertex & fragment shaders
// ========================================================

GLuint LoadShaders(const std::string &vertex_file_path,const std::string &fragment_file_path){

  std::cout << "load shaders" << std::endl;

  // Create the shaders
  GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
  GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
  
  // Read the Vertex Shader code from the file
  std::string VertexShaderCode;
  std::ifstream VertexShaderStream(vertex_file_path.c_str(), std::ios::in);
  if (VertexShaderStream.is_open()){
    std::string Line = "";
    while(getline(VertexShaderStream, Line))
      VertexShaderCode += "\n" + Line;
    VertexShaderStream.close();
  } else {
    std::cerr << "ERROR: cannot open " << vertex_file_path << std::endl;
    exit(0);
  }
  // Read the Fragment Shader code from the file
  std::string FragmentShaderCode;
  std::ifstream FragmentShaderStream(fragment_file_path.c_str(), std::ios::in);
  if(FragmentShaderStream.is_open()){
    std::string Line = "";
    while(getline(FragmentShaderStream, Line))
      FragmentShaderCode += "\n" + Line;
    FragmentShaderStream.close();
  } else {
    std::cerr << "ERROR: cannot open " << vertex_file_path << std::endl;
    exit(0);
  }
  
  GLint Result = GL_FALSE;
  int InfoLogLength;
  
  // Compile Vertex Shader
  std::cout << "Compiling shader : " << vertex_file_path << std::endl;
  char const * VertexSourcePointer = VertexShaderCode.c_str();
  glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
  glCompileShader(VertexShaderID);
  // Check Vertex Shader
  glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
  glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
  if ( InfoLogLength > 0 ){
    std::vector<char> VertexShaderErrorMessage(InfoLogLength+1);
    glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
    std::cerr << "VERTEX SHADER ERROR: " << std::string(VertexShaderErrorMessage.begin(),
                                                        VertexShaderErrorMessage.end()) << std::endl;
  }
  
  // Compile Fragment Shader
  std::cout << "Compiling shader : " << fragment_file_path << std::endl;
  char const * FragmentSourcePointer = FragmentShaderCode.c_str();
  glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
  glCompileShader(FragmentShaderID);
  // Check Fragment Shader
  glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
  glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
  if ( InfoLogLength > 0 ){
    std::vector<char> FragmentShaderErrorMessage(InfoLogLength+1);
    glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
    std::cerr << "FRAGMENT SHADER ERROR: " << std::string(FragmentShaderErrorMessage.begin(),
                                                          FragmentShaderErrorMessage.end()) << std::endl;
  }
  
  // Link the program
  std::cout << "Linking program" << std::endl;
  GLuint ProgramID = glCreateProgram();
  glAttachShader(ProgramID, VertexShaderID);
  glAttachShader(ProgramID, FragmentShaderID);
  glLinkProgram(ProgramID);
  // Check the program
  glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
  glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
  if ( InfoLogLength > 0 ){
    std::vector<char> ProgramErrorMessage(InfoLogLength+1);
    glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
    std::cerr << "SHADER PROGRAM ERROR: " << std::string(ProgramErrorMessage.begin(),
                                                         ProgramErrorMessage.end()) << std::endl;
  }
  
  glDeleteShader(VertexShaderID);
  glDeleteShader(FragmentShaderID);
  
  return ProgramID;
}

// ========================================================
// Functions related to error handling
// ========================================================

void GLCanvas::error_callback(int error, const char* description) {
  std::cerr << "ERROR CALLBACK: " << description << std::endl;
}

std::string WhichGLError(GLenum &error) {
  switch (error) {
  case GL_NO_ERROR:
    return "NO_ERROR";
  case GL_INVALID_ENUM:
    return "GL_INVALID_ENUM";
  case GL_INVALID_VALUE:
    return "GL_INVALID_VALUE";
  case GL_INVALID_OPERATION:
    return "GL_INVALID_OPERATION";
  case GL_INVALID_FRAMEBUFFER_OPERATION:
    return "GL_INVALID_FRAMEBUFFER_OPERATION";
  case GL_OUT_OF_MEMORY:
    return "GL_OUT_OF_MEMORY";
  case GL_STACK_UNDERFLOW:
    return "GL_STACK_UNDERFLOW";
  case GL_STACK_OVERFLOW:
    return "GL_STACK_OVERFLOW";
  default:
    return "OTHER GL ERROR";
  }
}

int HandleGLError(const std::string &message, bool ignore) {
  GLenum error;
  int i = 0;
  while ((error = glGetError()) != GL_NO_ERROR) {
    if (!ignore) {
      if (message != "") {
	std::cerr << "[" << message << "] ";
      }
      std::cerr << "GL ERROR(" << i << ") " << WhichGLError(error) << std::endl;
    }
    i++;
  }
  if (i == 0) return 1;
  return 0;
}

// ========================================================
// ========================================================


