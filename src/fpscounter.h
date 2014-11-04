#include "glCanvas.h"
#include "argparser.h"

// Includes for fps counter
#include <unistd.h>

#if _WIN32
#include <Winsock2.h>
#else
#include <sys/time.h>
#endif


// Call inside animate loop
// Will print out frames
void fpsCounter( const bool & active, float & timer ){

  // Requires you active == true if you want to keep track of time
  // timer will return updated time
  
#if _WIN32    
  static int last_tick_count;
  static int last_fps_count;
  static int frames = -1;
  if (frames == -1) {    
    last_fps_count = last_tick_count = GetTickCount();
    frames = 0;
  }
  if (active) {
    frames++;
    int this_tick_count = GetTickCount();
    // difference in milliseconds
    double diff = 0.001*(this_tick_count-last_tick_count);
    last_tick_count = this_tick_count;
    timer += diff;
    double diff_fps_time = 0.001*(this_tick_count-last_fps_count);
    if (diff_fps_time > 1.00) {      
      float fps = frames / float(diff_fps_time);
      std::cout << "fps: " << fps << std::endl;
      frames = 0;
      last_fps_count = this_tick_count;
    }
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
  if (active) {
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
    timer += diff;
  } else {
    gettimeofday(&last_time, NULL);
    last_fps_time = last_time;
    frames = 0;
    // if we aren't animating the light source, avoid busy-waiting!
    //usleep (100);
  }
#endif
  
}
