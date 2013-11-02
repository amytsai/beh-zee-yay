#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <Eigen/Dense>
#include <time.h>
#include <math.h>


#define PI 3.14159265  // Should be used from mathlib
inline float sqr(float x) { return x*x; }

using namespace std;

//****************************************************
// Some Classes
//****************************************************

class Viewport;

class Viewport {
  public:
    int w, h; // width and height
};

//****************************************************
// Global Variables
//****************************************************
Viewport	viewport;
float parameter;
bool adaptive = false;

//****************************************************
// Simple init function
//****************************************************
void initScene(){

  // Nothing to do here for this simple example.

}


//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
  viewport.w = w;
  viewport.h = h;

  glViewport (0,0,viewport.w,viewport.h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, viewport.w, 0, viewport.h);

}

//****************************************************
// A routine to set a pixel by drawing a GL point.  This is not a
// general purpose routine as it assumes a lot of stuff specific to
// this example.
//****************************************************

void setPixel(int x, int y, GLfloat r, GLfloat g, GLfloat b) {
  if (isToon) {
    glColor3f(toonify(r), toonify(g), toonify(b));
  } else {
    glColor3f(r, g, b);
  } 
  glVertex2f(x + 0.5, y + 0.5);   // The 0.5 is to target pixel
  // centers 
  // Note: Need to check for gap
  // bug on inst machines.
}

//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void myDisplay() {

  glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer

  glMatrixMode(GL_MODELVIEW);			        // indicate we are specifying camera transformations
  glLoadIdentity();				        // make sure transformation is "zero'd"


  // Start drawing
  if(isTor) {
	  torus(viewport.w / 2.0 , viewport.h / 2.0 , innerRad, outerRad);
  }
  else {
		circle(viewport.w / 2.0 , viewport.h / 2.0 , min(viewport.w, viewport.h) / 3.0);
  }
  glFlush();
  glutSwapBuffers();					// swap buffers (we earlier set double buffer)
}

//****************************************************
// Parse Scene File
//****************************************************

void loadScene(std::string file) {
  cout << "loading Scene .. "<< endl;
  ifstream inpfile(file.c_str());
  if(!inpfile.is_open()) {
    std::cout << "Unable to open file" << std::endl;
  } else {
    std::string line;

    while(inpfile.good()) {
      vector<string> splitline;
      string buf;

      getline(inpfile,line);
      stringstream ss(line);

      while (ss >> buf) {
        splitline.push_back(buf);
      }
      //Ignore blank lines
      if(splitline.size() == 0) {
        continue;
      }

      //Ignore comments
      if(splitline[0][0] == '#') {
        continue;
      }

      //size width height
      else if(!splitline[0].compare("size")) {
        width = atoi(splitline[1].c_str());
        height = atoi(splitline[2].c_str());
        bitmap = FreeImage_Allocate(width, height, 24);
        printf("Outputting image of size: %d x %d\n", width, height);
      }
      //maxdepth depth
      else if(!splitline[0].compare("maxdepth")) {
        maxdepth = atoi(splitline[1].c_str());
        printf("Raytracing with maxdepth = %d\n", maxdepth);
      }
      //output filename
      else if(!splitline[0].compare("output")) {
        filename = splitline[1];
        printf("Writing to file: %s\n", filename.c_str());
      } else {
        std::cerr << "Unknown command: " << splitline[0] << std::endl;
      }
    }
  }

}

//****************************************************
// MAIN
//****************************************************
int main(int argc, char *argv[]) {
  loadScene(argv[1]);
  if (argc < 3) {
    cout << "Not enough arguments" << endl;
    exit(EXIT_FAILURE);
  } else {
    parameter = argv[2];
  }

  if (argc == 3 ) {
      if (strcmp(argv[3], "-a") == 0) {
        adaptive = true;
      } else {
        cout << "Command line argument not found" << endl;
        exit(EXIT_FAILURE);
      }
  }

  //This initializes glut
  glutInit(&argc, argv);

  //This tells glut to use a double-buffered window with red, green, and blue channels 
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

  // Initalize theviewport size
  viewport.w = 400;
  viewport.h = 400;

  //The size and position of the window
  glutInitWindowSize(viewport.w, viewport.h);
  glutInitWindowPosition(0,0);
  glutCreateWindow(argv[0]);

  initScene();							// quick function to set up scene

  glutDisplayFunc(myDisplay);        // function to run when its time to draw something
  glutReshapeFunc(myReshape);        // function to run when the window gets resized

  glutMainLoop();							// infinite loop that will keep drawing and resizing
  // and whatever else

  return 0;
}








